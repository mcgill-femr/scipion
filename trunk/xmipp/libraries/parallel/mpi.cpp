/***************************************************************************
 * Authors:     J.M. de la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "mpi.h"

/** Function for a thread waiting on master MPI node distributing tasks.
 * This function will be called for one thread from the constructor
 * of the MpiTaskDistributor in the master.
 */
void __threadMpiMasterDistributor(ThreadArgument &arg)
{
    MpiTaskDistributor * distributor = (MpiTaskDistributor*) arg.workClass;
    int size = distributor->node->size;
    longint workBuffer[3];
    MPI_Status status;
    int finalizedWorkers = 0;

    while (finalizedWorkers < size - 1)
    {
        //wait for request form workers
        MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        workBuffer[0] = distributor->ThreadTaskDistributor::distribute(
                workBuffer[1], workBuffer[2]) ? 1 : 0;
        if (workBuffer[0] == 0) //no more jobs
            finalizedWorkers++;
        //send work
        MPI_Send(workBuffer, 3, MPI_LONG_LONG_INT, status.MPI_SOURCE, TAG_WORK,
                MPI_COMM_WORLD);
    }
}

MpiTaskDistributor::MpiTaskDistributor(longint nTasks, longint bSize,
        MpiNode *node) :
    ThreadTaskDistributor(nTasks, bSize)
{
    this->node = node;
    //if master create distribution thread
    if (node->isMaster())
    {
        manager = new ThreadManager(1, (void*) this);
        manager->runAsync(__threadMpiMasterDistributor);
    }
}

MpiTaskDistributor::~MpiTaskDistributor()
{
    if (node->isMaster())
    {
        manager->wait();
        delete manager;
    }
}

bool MpiTaskDistributor::distribute(longint &first, longint &last)
{
    if (node->isMaster())
        return ThreadTaskDistributor::distribute(first, last);

    //If not master comunicate with master thread
    //to get tasks
    //workBuffer[0] = 0 if no more jobs, 1 otherwise
    //workBuffer[1] = first
    //workBuffer[2] = last
    longint workBuffer[3];
    MPI_Status status;
    //any message from the master, is tag is TAG_STOP then stop
    MPI_Send(0, 0, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Recv(workBuffer, 3, MPI_LONG_LONG_INT, 0, TAG_WORK, MPI_COMM_WORLD,
            &status);

    first = workBuffer[1];
    last = workBuffer[2];

    return (workBuffer[0] == 1);
}

FileTaskDistributor::FileTaskDistributor(longint nTasks, longint bSize,
        MpiNode * node) :
    ThreadTaskDistributor(nTasks, bSize)
{
    fileCreator = false;

    if (node == NULL || node->isMaster())
    {
        fileCreator = true;
        //not using mpi or in master node
        fileCreator = true;
        strcpy(lockFilename, "pijol_XXXXXX");
        if ((lockFile = mkstemp(lockFilename)) == -1)
        {
            perror("FileTaskDistributor::Error generating tmp lock file");
            exit(1);
        }
        else
            close(lockFile);
        createLockFile();
    }
    //if using mpi broadcast the filename from master to slaves
    if (node != NULL)
        MPI_Bcast(lockFilename, L_tmpnam, MPI_CHAR, 0, MPI_COMM_WORLD);
    loadLockFile();
}

FileTaskDistributor::~FileTaskDistributor()
{
    close(lockFile);
    if (fileCreator && remove(lockFilename) == -1)
        perror("FileTaskDistributor: error deleting lock file");
}

void FileTaskDistributor::createLockFile()
{

    int buffer[] = { numberOfTasks, assignedTasks, blockSize };

    if ((lockFile = open(lockFilename, O_CREAT | O_RDWR | O_TRUNC, S_IRUSR
            | S_IWUSR | S_IRGRP | S_IROTH)) == -1)
    {
        perror("FileTaskDistributor::createLockFile: Error opening lock file");
        exit(1);
    }

    if (write(lockFile, buffer, 3 * sizeof(int)) == -1)
    {
        perror(
                "FileTaskDistributor::createLockFile: Error writing to lock file");
        exit(1);
    }

    writeVars();
}//function createLockFile

void FileTaskDistributor::loadLockFile()
{
    if ((lockFile = open(lockFilename, O_RDWR)) == -1)
    {
        perror("FileTaskDistributor::loadLockFile: Error opening lock file");
        exit(1);
    }
    readVars();
}

void FileTaskDistributor::readVars()
{
    lseek(lockFile, 0, SEEK_SET);
    read(lockFile, &numberOfTasks, sizeof(long long int));
    read(lockFile, &assignedTasks, sizeof(long long int));
    read(lockFile, &blockSize, sizeof(long long int));
}

void FileTaskDistributor::writeVars()
{
    lseek(lockFile, 0, SEEK_SET);
    write(lockFile, &numberOfTasks, sizeof(long long int));
    write(lockFile, &assignedTasks, sizeof(long long int));
    write(lockFile, &blockSize, sizeof(long long int));
}

void FileTaskDistributor::lock()
{
    ThreadTaskDistributor::lock();
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_LOCK, 0);
    readVars();
}

void FileTaskDistributor::unlock()
{
    writeVars();
    lseek(lockFile, 0, SEEK_SET);
    lockf(lockFile, F_ULOCK, 0);
    ThreadTaskDistributor::unlock();
}

//------------ MPI ---------------------------
MpiNode::MpiNode(int &argc, char ** argv)
{
    //MPI Initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

MpiNode::~MpiNode()
{
    MPI_Finalize();
}

bool MpiNode::isMaster() const
{
    return rank == 0;
}

void MpiNode::barrierWait()
{
  MPI_Barrier(MPI_COMM_WORLD);
}
