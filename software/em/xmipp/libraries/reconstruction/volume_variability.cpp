/***************************************************************************
 *
 * Authors:     Javier Vargas (javier.vargasbalbuena@mcgill.ca)
 *
 *
 * Department of Anatomy and Cell Biology, McGill University
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
 *  e-mail address 'javier.vargasbalbuena@mcgill.ca'
 ***************************************************************************/


/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *              Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Jose Roman Bilbao-Castro (jrbcast@ace.ual.es)
 *              Vahid Abrishami (vabrishami@cnb.csic.es)
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


#include "volume_variability.h"
#include <iostream>

// Define params
void ProgVolVariability::defineParams()
{
	//usage
    addUsageLine("Generate variability maps using direct Fourier interpolation with arbitrary geometry ");
    addUsageLine("from projections and a previous 3D reconstruction. Kaisser-windows are used for interpolation in Fourier space.");
    //params
    addParamsLine("   -i <md_file>                : Metadata file with input projections");
    addParamsLine("   --vol <volume_file=\"rec_fourier.vol\"> : Filename for input volume");
    addParamsLine("  [-o <volume_file=\"rec_fourier.vol\">]   : Filename for output variability volume");
	addParamsLine("  [--mask <vol_file=\"\">]     : Binary Mask defining the macromolecule");
    addParamsLine("  [--iter <iterations=1>]      : Number of iterations for weight correction");
    addParamsLine("  [--sym <symfile=c1>]              : Enforce symmetry in projections");
    addParamsLine("  [--padding <proj=2.0> <vol=2.0>]  : Padding used for projections and volume");
    addParamsLine("  [--prepare_fsc <fscfile>]      : Filename root for FSC files");
    addParamsLine("  [--min_resolution <p=0.0>]     	: Min resolution (default=0.0)");
    addParamsLine("  [--max_resolution <p=0.5>]     : Max resolution (Nyquist=0.5)");
    addParamsLine("  [--weight]                     : Use weights stored in the image metadata");
    addParamsLine("  [--thr <threads=1> <rows=1>]   : Number of concurrent threads and rows processed at time by a thread");
    addParamsLine("  [--thrFFT <threadsFFT=1>]   	: Number of concurrent threads for FFT calculations");
    addParamsLine("  [--blob <radius=1.9> <order=0> <alpha=15>] : Blob parameters");
    addParamsLine("                                 : radius in pixels, order of Bessel function in blob and parameter alpha");
    addParamsLine("  [--useCTF]                     : Use CTF information if present");
    addParamsLine("  [--sampling <Ts=1>]            : sampling rate of the input images in Angstroms/pixel");
    addParamsLine("                                 : It is only used when correcting for the CTF");
    addParamsLine("  [--phaseFlipped]               : Give this flag if images have been already phase flipped");
    addParamsLine("  [--minCTF <ctf=0.01>]          : Minimum value of the CTF that will be inverted");
    addParamsLine("                                 : CTF values (in absolute value) below this one will not be corrected");
    addParamsLine("  [--iterMC <iterations=100>]    : Number of Monte Carlo iterations");

 //   addExampleLine("For reconstruct enforcing i3 symmetry and using stored weights:", false);
 //   addExampleLine("   xmipp_reconstruct_fourier  -i reconstruction.sel --sym i3 --weight");
 //   JV QUESTIONS;
    // NiterWeight por defecto esta a 1, pero en muchos puntos se hace la comparacion con NiterWeight == 0, y me da la sensacion
    // 			   que NiterWeight == 0, deberia ser NiterWeight == 1 mirar por ejem lineas 186-189
    // Tiene sentido usar numThreads para la FFT y forzar que numThreads es 1 para la parte de reconstruccion?
    // usar prepare_fsc para calcular FSC con los halves ids bien. Usar metadata _rln_halfid

    //Cosas probar: no distorsional Vin() por la FT del blob ? mejora algo?



    //Tests: las trasnformaciones de Vin y fftIn en ProduceSideInfo estan bien? Comprobarlo quitando toda average y ver que sale bien.
}

// Read arguments ==========================================================
void ProgVolVariability::readParams()
{
    fn_sel = getParam("-i");
    fn_out = getParam("-o");
    fn_sym = getParam("--sym");

    //~JV
    fn_vol = getParam("--vol");
    fn_mask = getParam("--mask");
    //~JV

    if(checkParam("--prepare_fsc"))
        fn_fsc = getParam("--prepare_fsc");
    do_weights = checkParam("--weight");
    padding_factor_proj = getDoubleParam("--padding", 0);
    padding_factor_vol = getDoubleParam("--padding", 1);
    blob.radius   = getDoubleParam("--blob", 0);
    blob.order    = getIntParam("--blob", 1);
    blob.alpha    = getDoubleParam("--blob", 2);
    minResolution = getDoubleParam("--min_resolution");
    maxResolution = getDoubleParam("--max_resolution");
    numThreads = getIntParam("--thr");
    thrWidth = getIntParam("--thr", 1);
    numThreadsFFT = getIntParam("--thrFFT");
    NiterWeight = getIntParam("--iter");
    NiterMC = getIntParam("--iterMC");
    useCTF = checkParam("--useCTF");
    phaseFlipped = checkParam("--phaseFlipped");
    minCTF = getDoubleParam("--minCTF");
    if (useCTF)
        Ts=getDoubleParam("--sampling");

}

// Show ====================================================================
void ProgVolVariability::show()
{
    if (verbose > 0)
    {
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Variability estimation method using Kaiser windows as interpolators"   << std::endl;
        std::cout << " =====================================================================" << std::endl;
        std::cout << " Input selfile             		: "  << fn_sel << std::endl;
        std::cout << " Input reconstructed average map  : "  << fn_vol << std::endl;
        std::cout << " Input mask map  					: "  << fn_mask << std::endl;
        std::cout << " padding_factor_proj       		: "  << padding_factor_proj << std::endl;
        std::cout << " padding_factor_vol        		: "  << padding_factor_vol << std::endl;
        std::cout << " Output variability map      		: "  << fn_out << std::endl;
        if (fn_sym != "")
            std::cout << " Symmetry file for projections : "  << fn_sym << std::endl;
        if (fn_fsc != "")
            std::cout << " File root for FSC files: " << fn_fsc << std::endl;
        if (do_weights)
            std::cout << " Use weights stored in the image headers or doc file" << std::endl;
        else
            std::cout << " Do NOT use weights" << std::endl;
        if (useCTF)
            std::cout << "Using CTF information" << std::endl
            << "Sampling rate: " << Ts << std::endl
            << "Phase flipped: " << phaseFlipped << std::endl
            << "Minimum CTF: " << minCTF << std::endl;
        std::cout << "\n Interpolation Function"
        << "\n   blrad                 : "  << blob.radius
        << "\n   blord                 : "  << blob.order
        << "\n   blalpha               : "  << blob.alpha
        << "\n sampling_rate           : "  << Ts
        << "\n min_resolution          : "  << minResolution
        << "\n max_resolution          : "  << maxResolution
        << "\n -----------------------------------------------------------------" << std::endl;

    }
}

// Main routine ------------------------------------------------------------
void ProgVolVariability::run()
{
    show();
    produceSideinfo();
    // Process all images in the selfile
    if (verbose)
    {
        if (NiterWeight!=0)
            init_progress_bar(NiterWeight*SF.size());
        else
            init_progress_bar(SF.size());
    }
    // Create threads stuff
    barrier_init( &barrier, numThreads+1 );
    pthread_mutex_init( &workLoadMutex, NULL );
    statusArray = NULL;
    th_ids = (pthread_t *)malloc( numThreads * sizeof( pthread_t));
    th_args = (ImageThreadParams *) malloc ( numThreads * sizeof( ImageThreadParams ) );

    // Create threads
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
    {
        // Passing parameters to each thread
        th_args[nt].parent = this;
        th_args[nt].myThreadID = nt;
        th_args[nt].selFile = new MetaData(SF);
        pthread_create( (th_ids+nt) , NULL, processImageThread, (void *)(th_args+nt) );
    }

    //Computing interpolated volume

    processImages(0, SF.size() - 1, !fn_fsc.empty(), false);

    // Correcting the weights
    correctWeight();

    //Saving the volume
    finishComputations(fn_out);

    threadOpCode = EXIT_THREAD;

    // Waiting for threads to finish
    barrier_wait( &barrier );
    for ( int nt = 0 ; nt < numThreads ; nt ++ )
        pthread_join(*(th_ids+nt), NULL);
    barrier_destroy( &barrier );

    // Deallocate resources.
    if ( statusArray != NULL)
    {
        free(statusArray);
    }
    for (int nt=0; nt<numThreads ;nt++)
    {
    	delete(th_args[nt].selFile);
    }
    free(th_ids);
    free(th_args);
}


void ProgVolVariability::produceSideinfo()
{
    // Translate the maximum resolution to digital frequency
    // maxResolution=sampling_rate/maxResolution;
    maxResolution2=maxResolution*maxResolution;
    minResolution2=minResolution*minResolution;

    // Read the input images
    SF.read(fn_sel);
    SF.removeDisabled();

    // Ask for memory for the output volume and its Fourier transform
    size_t objId = SF.firstObject();
    FileName fnImg;
    SF.getValue(MDL_IMAGE,fnImg,objId);
    Image<double> I;
    I.read(fnImg, HEADER);
    int Ydim=YSIZE(I());
    int Xdim=XSIZE(I());
    if (Ydim!=Xdim)
        REPORT_ERROR(ERR_MULTIDIM_SIZE,"This algorithm only works for squared images");
    imgSize=Xdim;

    volPadSizeX = volPadSizeY = volPadSizeZ=(int)(Xdim*padding_factor_vol);

    //use threads for volume inverse fourier transform, plan is created in setReal()
    transformerVol.setThreadsNumber(numThreadsFFT);

    Vout().initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    transformerVol.cleanup();
    transformerVol.setReal(Vout());

#ifdef DEBUG
	std::cout << A3D_ELEM(fftVin, 0, 0, 0) << std::endl;
	std::cout << A3D_ELEM(fftVin, 1, 1, 1) << std::endl;
	std::cout << A3D_ELEM(fftVin, 2, 1, 1) << std::endl;
#endif

    Vout().clear(); // Free the memory so that it is available for FourierWeights
    transformerVol.getFourierAlias(VoutFourier);
    VoutFourier.initZeros();
    FourierWeights.initZeros(VoutFourier);

    // Ask for memory for the padded images
    size_t paddedImgSize=(size_t)(Xdim*padding_factor_proj);
    paddedImg.resize(paddedImgSize,paddedImgSize);
    paddedImg.setXmippOrigin();
    transformerImg.setReal(paddedImg);

    // Build a table of blob values
    blobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    fourierBlobTableSqrt.resize(BLOB_TABLE_SIZE_SQRT);
    Fourier_blob_table.resize(BLOB_TABLE_SIZE_SQRT);

    struct blobtype blobFourier,blobnormalized;
    blobFourier=blob;
    //Sjors 18aug10 blobFourier.radius/=(padding_factor_proj*Xdim);
    blobFourier.radius/=(padding_factor_vol*Xdim);
    blobnormalized=blob;
    blobnormalized.radius/=((double)padding_factor_proj/padding_factor_vol);
    double deltaSqrt     = (blob.radius*blob.radius) /(BLOB_TABLE_SIZE_SQRT-1);
    double deltaFourier  = (sqrt(3.)*Xdim/2.)/(BLOB_TABLE_SIZE_SQRT-1);

    // The interpolation kernel must integrate to 1
    double iw0 = 1.0 / blob_Fourier_val(0.0, blobnormalized);
    //Sjors 18aug10 double padXdim3 = padding_factor_proj * Xdim;
    double padXdim3 = padding_factor_vol * Xdim;
    padXdim3 = padXdim3 * padXdim3 * padXdim3;
    double blobTableSize = blob.radius*sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    //***
    //Following commented line seems to be the right thing but I do not understand it
    //double fourierBlobTableSize = (sqrt(3.)*Xdim*Xdim/2.)*blobFourier.radius *sqrt(1./ (BLOB_TABLE_SIZE_SQRT-1));
    FOR_ALL_ELEMENTS_IN_MATRIX1D(blobTableSqrt)

    {
        //use a r*r sample instead of r
        //DIRECT_VEC_ELEM(blob_table,i)         = blob_val(delta*i, blob)  *iw0;
        VEC_ELEM(blobTableSqrt,i)    = blob_val(blobTableSize*sqrt((double)i), blob)  *iw0;
        //***
        //DIRECT_VEC_ELEM(fourierBlobTableSqrt,i) =
        //     blob_Fourier_val(fourierBlobTableSize*sqrt(i), blobFourier)*padXdim3  *iw0;
        VEC_ELEM(Fourier_blob_table,i) =
            blob_Fourier_val(deltaFourier*i, blobFourier)*padXdim3  *iw0;
        //#define DEBUG
#ifdef DEBUG

        std::cout << VEC_ELEM(Fourier_blob_table,i)
        << " " << VEC_ELEM(fourierBlobTableSqrt,i)
        << std::endl;
#endif
  #undef DEBUG

    }
    //iDelta        = 1/delta;
    iDeltaSqrt    = 1/deltaSqrt;
    iDeltaFourier = 1/deltaFourier;

    // Get symmetries
    Matrix2D<double>  Identity(3,3);
    Identity.initIdentity();
    R_repository.push_back(Identity);
    if (fn_sym != "")
    {
        SymList SL;
        SL.readSymmetryFile(fn_sym);
        for (int isym = 0; isym < SL.symsNo(); isym++)
        {
            Matrix2D<double>  L(4, 4), R(4, 4);
            SL.getMatrices(isym, L, R);
            R.resize(3, 3);
            R_repository.push_back(R);
        }
    }

    // ~JV
    //Read the input average map
    //Input volume Vin
    Image<double> Vin;
    Vin.read(fn_vol);
    Vin().setXmippOrigin();

    //JV: Distort Vin() by the Fourier Transform of the Blob (see lines starting 1356 in finishComputations)
    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
    pad_relation = (pad_relation * pad_relation * pad_relation);

    MultidimArray<double> &mVin=Vin();

    double ipad_relation=1.0/pad_relation;
    double meanFactor2=0;

    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
    {
        double radius=sqrt((double)(k*k+i*i+j*j));
        double aux=radius*iDeltaFourier;
        double factor = Fourier_blob_table(ROUND(aux));
        double factor2=(pow(Sinc(radius/(2*(imgSize))),2));
        if (NiterWeight!=0)
        {
            A3D_ELEM(mVin,k,i,j) *= (ipad_relation*factor2*factor);  //JV changed!
            meanFactor2+=factor2;
        }
        else
            A3D_ELEM(mVin,k,i,j) *= (ipad_relation*factor);
    }

    if (NiterWeight!=0)
    {
        meanFactor2/=MULTIDIM_SIZE(mVin);
        FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
        A3D_ELEM(mVin,k,i,j) /= meanFactor2;  //JV changed!
    }
    //JV: Finish here Vin distortion by FT of Blob

    //Padding of Vin
    MultidimArray<double> localPaddedVin;
    localPaddedVin.initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    localPaddedVin.setXmippOrigin();

    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVin)
    A3D_ELEM(localPaddedVin,k,i,j) = A3D_ELEM(mVin,k,i,j);


#ifdef DEBUG
    Image<double> save;
    save().alias( localPaddedVin );
    save.write((std::string)"Vin.vol");
#endif

    CenterFFT(localPaddedVin,true);
    FourierTransformer transformerVol2;
    transformerVol2.setThreadsNumber(numThreadsFFT);
    fftVin.initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);
    transformerVol2.FourierTransform(localPaddedVin, fftVin, true);
    transformerVol2.clear();

    //JV See lines starting 596, This is because the normalizing factors in the Fourier Transform?
    double corr3D_2D=(imgSize* pow(padding_factor_vol,3.)) / pow(padding_factor_proj,2.);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(fftVin)
    A3D_ELEM(fftVin,k,i,j)*=corr3D_2D;

    randomize_random_generator();
    localPaddedVin.clear();

#ifdef DEBUG
    Image< std::complex<double> > save;
    save().alias(fftVin);
    save.write((std::string)"Vin.vol");
#endif

    // ~JV

}

void * ProgVolVariability::processImageThread( void * threadArgs )
{

    ImageThreadParams * threadParams = (ImageThreadParams *) threadArgs;
    ProgVolVariability * parent = threadParams->parent;
    barrier_t * barrier = &(parent->barrier);

    int minSeparation;

    if ( (int)ceil(parent->blob.radius) > parent->thrWidth )
        minSeparation = (int)ceil(parent->blob.radius);
    else
        minSeparation = parent->thrWidth;

    minSeparation+=1;

    Matrix2D<double>  localA(3, 3), localAinv(3, 3);
    MultidimArray< std::complex<double> > localPaddedFourier;
    MultidimArray<double> localPaddedImg;
    FourierTransformer localTransformerImg;

    std::vector<size_t> objId;

    threadParams->selFile->findObjects(objId);
    ApplyGeoParams params;
    params.only_apply_shifts = true;
    MultidimArray<int> zWrapped(3*parent->volPadSizeZ),yWrapped(3*parent->volPadSizeY),xWrapped(3*parent->volPadSizeX),
    zNegWrapped, yNegWrapped, xNegWrapped;
    zWrapped.initConstant(-1);
    yWrapped.initConstant(-1);
    xWrapped.initConstant(-1);
    zWrapped.setXmippOrigin();
    yWrapped.setXmippOrigin();
    xWrapped.setXmippOrigin();
    zNegWrapped=zWrapped;
    yNegWrapped=yWrapped;
    xNegWrapped=xWrapped;

    MultidimArray<double> x2precalculated(XSIZE(xWrapped)), y2precalculated(XSIZE(yWrapped)), z2precalculated(XSIZE(zWrapped));
    x2precalculated.initConstant(-1);
    y2precalculated.initConstant(-1);
    z2precalculated.initConstant(-1);
    x2precalculated.setXmippOrigin();
    y2precalculated.setXmippOrigin();
    z2precalculated.setXmippOrigin();

    bool hasCTF=(threadParams->selFile->containsLabel(MDL_CTF_MODEL) || threadParams->selFile->containsLabel(MDL_CTF_DEFOCUSU)) &&
                parent->useCTF;
    if (hasCTF)
    {
        threadParams->ctf.enable_CTF=true;
        threadParams->ctf.enable_CTFnoise=false;
    }
    do
    {
        barrier_wait( barrier );

        switch ( parent->threadOpCode )
        {
        case PRELOAD_IMAGE:
            {
                threadParams->read = 0;

                if ( threadParams->imageIndex >= 0 )
                {
                    // Read input image
                    double rot, tilt, psi, weight;
                    Projection proj;

                    //Read projection from selfile, read also angles and shifts if present
                    //but only apply shifts
                    proj.readApplyGeo(*(threadParams->selFile), objId[threadParams->imageIndex], params);
                    rot  = proj.rot();
                    tilt = proj.tilt();
                    psi  = proj.psi();
                    weight = proj.weight();
                    if (hasCTF)
                    {
                        threadParams->ctf.readFromMetadataRow(*(threadParams->selFile),objId[threadParams->imageIndex]);
                        // threadParams->ctf.Tm=threadParams->parent->Ts;
                        threadParams->ctf.produceSideInfo();
                    }

                    threadParams->weight = 1.;

                    if(parent->do_weights)
                        threadParams->weight = weight;
                    else if (!parent->do_weights)
                    {
                        weight=1.0;
                    }
                    else if (weight==0.0) // JV I think that the flux never enters here!
                    {
                        threadParams->read = 2;
                        break;
                    }

                    // Copy the projection to the center of the padded image
                    // and compute its Fourier transform
                    proj().setXmippOrigin();
                    size_t localPaddedImgSize=(size_t)(parent->imgSize*parent->padding_factor_proj);
                    if (threadParams->reprocessFlag) //JV if I use this flag to recalculate the 3D reconstruction this line should be removed
                        localPaddedFourier.initZeros(localPaddedImgSize,localPaddedImgSize/2+1);
                    else
                    {
                        localPaddedImg.initZeros(localPaddedImgSize,localPaddedImgSize);
                        localPaddedImg.setXmippOrigin();
                        const MultidimArray<double> &mProj=proj();
                        FOR_ALL_ELEMENTS_IN_ARRAY2D(mProj)
                        A2D_ELEM(localPaddedImg,i,j)=A2D_ELEM(mProj,i,j);
                        // COSS A2D_ELEM(localPaddedImg,i,j)=weight*A2D_ELEM(mProj,i,j);
                        CenterFFT(localPaddedImg,true);

                        // Fourier transformer for the images
                        localTransformerImg.setReal(localPaddedImg);
                        localTransformerImg.FourierTransform();
                        localTransformerImg.getFourierAlias(localPaddedFourier);
                    }

                    // Compute the coordinate axes associated to this image
                    Euler_angles2matrix(rot, tilt, psi, localA);
                    localAinv=localA.transpose();

                    threadParams->localweight = weight;
                    threadParams->localAInv = &localAinv;
                    threadParams->localPaddedFourier = &localPaddedFourier;
                    //#define DEBUG22
#ifdef DEBUG22

                    {//CORRECTO

                        if(threadParams->myThreadID%1==0)
                        {
                            proj.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                       integerToString(threadParams->imageIndex) + "proj.spi");

                            ImageXmipp save44;
                            save44()=localPaddedImg;
                            save44.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_img.spi");

                            FourierImage save33;
                            save33()=localPaddedFourier;
                            save33.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "local_padded_fourier.spi");
                            FourierImage save22;
                            //save22()=*paddedFourier;
                            save22().alias(*(threadParams->localPaddedFourier));
                            save22.write((std::string) integerToString(threadParams->myThreadID)  + "_" +\
                                         integerToString(threadParams->imageIndex) + "_padded_fourier.spi");
                        }

                    }
#endif
                    #undef DEBUG22

                    threadParams->read = 1;
                }
                break;
            }

        case EXIT_THREAD:
            return NULL;

        case PROCESS_WEIGHTS:
            {
            	//JV: Here the differences are that go comes here is the variance which is the square of std. We correct corr2D_3D in finishComputations better

                // Get a first approximation of the reconstruction
                double corr2D_3D=pow(parent->padding_factor_proj,2.)/
                                 (parent->imgSize* pow(parent->padding_factor_vol,3.));

                // Divide by Zdim because of the
                // the extra dimension added
                // and padding differences

                //#ifdef DEBUG // without computing the variability the following should be the same
                //#endif

                MultidimArray<double> &mFourierWeights=parent->FourierWeights;
                for (int k=threadParams->myThreadID; k<=FINISHINGZ(mFourierWeights); k+=parent->numThreads)
                    for (int i=STARTINGY(mFourierWeights); i<=FINISHINGY(mFourierWeights); i++)
                        for (int j=STARTINGX(mFourierWeights); j<=FINISHINGX(mFourierWeights); j++)
                        {
                        	if (parent->NiterWeight==0)
//JV
//                        		A3D_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D;
                        		;
//JV
                        	else
                        	{
                        		double weight_kij=A3D_ELEM(mFourierWeights,k,i,j);
                        		if (1.0/weight_kij>ACCURACY)
//JV
//                        			A3D_ELEM(parent->VoutFourier,k,i,j)*=corr2D_3D*A3D_ELEM(mFourierWeights,k,i,j);
                        			A3D_ELEM(parent->VoutFourier,k,i,j)*=A3D_ELEM(mFourierWeights,k,i,j);
//JV
                        		else
                        			A3D_ELEM(parent->VoutFourier,k,i,j)=0;

                        	}

                        }
                break;
            }
        case PROCESS_IMAGE:
            {
                MultidimArray< std::complex<double> > *paddedFourier = threadParams->paddedFourier;
                if (threadParams->weight==0.0)
                    break;
                bool reprocessFlag = threadParams->reprocessFlag;
                int * statusArray = parent->statusArray;

                int minAssignedRow;
                int maxAssignedRow;
                bool breakCase;
                bool assigned;

                // Get the inverse of the sampling rate
                // double iTs=parent->padding_factor_proj/parent->Ts;
                double iTs=1.0/parent->Ts; // The padding factor is not considered here, but later when the indexes
                //                         // are converted to digital frequencies
                do
                {
                    minAssignedRow = -1;
                    maxAssignedRow = -1;
                    breakCase = false;
                    assigned = false;

                    do
                    {
                        pthread_mutex_lock( &(parent->workLoadMutex) );

                        if ( parent->rowsProcessed == YSIZE(*paddedFourier) )
                        {
                            pthread_mutex_unlock( &(parent->workLoadMutex) );
                            breakCase = true;
                            break;
                        }

                        for (size_t w = 0 ; w < YSIZE(*paddedFourier) ; w++ )
                        {
                            if ( statusArray[w]==0 )
                            {
                                assigned = true;
                                minAssignedRow = w;
                                maxAssignedRow = w+minSeparation-1;

                                if ( maxAssignedRow > (int)(YSIZE(*paddedFourier)-1) )
                                    maxAssignedRow = YSIZE(*paddedFourier)-1;

                                for ( int in = (minAssignedRow - minSeparation) ; in < (int)minAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }

                                for ( int in = minAssignedRow ; in <= (int)maxAssignedRow ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] == 0 )
                                        {
                                            statusArray[in] = -1;
                                            parent->rowsProcessed++;
                                        }
                                    }
                                }

                                for ( int in = maxAssignedRow+1 ; in <= (maxAssignedRow+minSeparation) ; in ++ )
                                {
                                    if ( ( in >= 0 ) && ( in < (int)YSIZE(*paddedFourier) ))
                                    {
                                        if ( statusArray[in] > -1 )
                                            statusArray[in]++;
                                    }
                                }

                                break;
                            }
                        }

                        pthread_mutex_unlock( &(parent->workLoadMutex) );
                    }
                    while ( !assigned );

                    if ( breakCase == true )
                    {
                        break;
                    }

                    Matrix2D<double> * A_SL = threadParams->symmetry;

                    // Loop over all Fourier coefficients in the padded image
                    Matrix1D<double> freq(3), gcurrent(3), real_position(3), contFreq(3);
                    Matrix1D<int> corner1(3), corner2(3);

                    // Some alias and calculations moved from heavy loops
                    double wCTF=1, wModulator=1.0;
                    double blobRadiusSquared = parent->blob.radius * parent->blob.radius;
                    double iDeltaSqrt = parent->iDeltaSqrt;
                    Matrix1D<double> & blobTableSqrt = parent->blobTableSqrt;
                    int xsize_1 = XSIZE(parent->VoutFourier) - 1;
                    int zsize_1 = ZSIZE(parent->VoutFourier) - 1;
                    MultidimArray< std::complex<double> > &VoutFourier=parent->VoutFourier;
                    MultidimArray<double> &fourierWeights = parent->FourierWeights;
                    //                    MultidimArray<double> &prefourierWeights = parent->preFourierWeights;
                    // Get i value for the thread
                    for (int i = minAssignedRow; i <= maxAssignedRow ; i ++ )
                    {
                        // Discarded rows can be between minAssignedRow and maxAssignedRow, check
                        if ( statusArray[i] == -1 )
                            for (int j=STARTINGX(*paddedFourier); j<=FINISHINGX(*paddedFourier); j++)
                            {
                                // Compute the frequency of this coefficient in the
                                // universal coordinate system
                                FFT_IDX2DIGFREQ(j,XSIZE(parent->paddedImg),XX(freq));
                                FFT_IDX2DIGFREQ(i,YSIZE(parent->paddedImg),YY(freq));
                                ZZ(freq)=0;
                                if (XX(freq)*XX(freq)+YY(freq)*YY(freq)>parent->maxResolution2)
                                    continue;
                                if (XX(freq)*XX(freq)+YY(freq)*YY(freq)<parent->minResolution2)
                                    continue;
                                wModulator=1.0;
                                if (hasCTF && !reprocessFlag)
                                {
                                    XX(contFreq)=XX(freq)*iTs;
                                    YY(contFreq)=YY(freq)*iTs;
                                    threadParams->ctf.precomputeValues(XX(contFreq),YY(contFreq));
                                    //wCTF=threadParams->ctf.getValueAt();
//JV
//                                    wCTF=threadParams->ctf.getValuePureNoKAt();
                                    wCTF=threadParams->ctf.getValuePureWithoutDampingAt();
//JV
                                    if (std::isnan(wCTF))
                                    {
                                    	if (i==0 && j==0)
                                    		wModulator=wCTF=1.0;
                                    	else
                                    		wModulator=wCTF=0.0;
                                    }
                                    if (fabs(wCTF)<parent->minCTF)
                                    {
                                        wModulator=fabs(wCTF);
                                        wCTF=SGN(wCTF);
                                    }

                                    else

                                    	wCTF=1.0/wCTF;

                                    if (parent->phaseFlipped)
                                        wCTF=fabs(wCTF);
                                }

                                SPEED_UP_temps012;
                                M3x3_BY_V3x1(freq,*A_SL,freq);

                                // Look for the corresponding index in the volume Fourier transform
                                DIGFREQ2FFT_IDX_DOUBLE(XX(freq),parent->volPadSizeX,XX(real_position));
                                DIGFREQ2FFT_IDX_DOUBLE(YY(freq),parent->volPadSizeY,YY(real_position));
                                DIGFREQ2FFT_IDX_DOUBLE(ZZ(freq),parent->volPadSizeZ,ZZ(real_position));

                                // Put a box around that coefficient
                                XX(corner1)=CEIL (XX(real_position)-parent->blob.radius);
                                YY(corner1)=CEIL (YY(real_position)-parent->blob.radius);
                                ZZ(corner1)=CEIL (ZZ(real_position)-parent->blob.radius);
                                XX(corner2)=FLOOR(XX(real_position)+parent->blob.radius);
                                YY(corner2)=FLOOR(YY(real_position)+parent->blob.radius);
                                ZZ(corner2)=FLOOR(ZZ(real_position)+parent->blob.radius);

#ifdef DEBUG

                                std::cout << "Idx Img=(0," << i << "," << j << ") -> Freq Img=("
                                << freq.transpose() << ") ->\n    Idx Vol=("
                                << real_position.transpose() << ")\n"
                                << "   Corner1=" << corner1.transpose() << std::endl
                                << "   Corner2=" << corner2.transpose() << std::endl;
#endif
                                // Loop within the box
                                double *ptrIn=(double *)&(A2D_ELEM(*paddedFourier, i,j));

                                // Some precalculations
                                for (int intz = ZZ(corner1); intz <= ZZ(corner2); ++intz)
                                {
                                    double z = intz - ZZ(real_position);
                                    A1D_ELEM(z2precalculated,intz)=z*z;
                                    if (A1D_ELEM(zWrapped,intz)<0)
                                    {
                                        int iz, izneg;
                                        fastIntWRAP(iz, intz, 0, zsize_1);
                                        A1D_ELEM(zWrapped,intz)=iz;
                                        int miz=-iz;
                                        fastIntWRAP(izneg, miz,0,zsize_1);
                                        A1D_ELEM(zNegWrapped,intz)=izneg;
                                    }
                                }
                                for (int inty = YY(corner1); inty <= YY(corner2); ++inty)
                                {
                                    double y = inty - YY(real_position);
                                    A1D_ELEM(y2precalculated,inty)=y*y;
                                    if (A1D_ELEM(yWrapped,inty)<0)
                                    {
                                        int iy, iyneg;
                                        fastIntWRAP(iy, inty, 0, zsize_1);
                                        A1D_ELEM(yWrapped,inty)=iy;
                                        int miy=-iy;
                                        fastIntWRAP(iyneg, miy,0,zsize_1);
                                        A1D_ELEM(yNegWrapped,inty)=iyneg;
                                    }
                                }
                                for (int intx = XX(corner1); intx <= XX(corner2); ++intx)
                                {
                                    double x = intx - XX(real_position);
                                    A1D_ELEM(x2precalculated,intx)=x*x;
                                    if (A1D_ELEM(xWrapped,intx)<0)
                                    {
                                        int ix, ixneg;
                                        fastIntWRAP(ix, intx, 0, zsize_1);
                                        A1D_ELEM(xWrapped,intx)=ix;
                                        int mix=-ix;
                                        fastIntWRAP(ixneg, mix,0,zsize_1);
                                        A1D_ELEM(xNegWrapped,intx)=ixneg;
                                    }
                                }

                                // Actually compute
                                for (int intz = ZZ(corner1); intz <= ZZ(corner2); ++intz)
                                {
                                    double z2 = A1D_ELEM(z2precalculated,intz);
                                    int iz=A1D_ELEM(zWrapped,intz);
                                    int izneg=A1D_ELEM(zNegWrapped,intz);

                                    for (int inty = YY(corner1); inty <= YY(corner2); ++inty)
                                    {
                                        double y2z2 = A1D_ELEM(y2precalculated,inty) + z2;
                                        if (y2z2 > blobRadiusSquared)
                                            continue;
                                        int iy=A1D_ELEM(yWrapped,inty);
                                        int iyneg=A1D_ELEM(yNegWrapped,inty);

                                        int	size1=YXSIZE(VoutFourier)*(izneg)+((iyneg)*XSIZE(VoutFourier));
                                        int	size2=YXSIZE(VoutFourier)*(iz)+((iy)*XSIZE(VoutFourier));
                                        int	fixSize=0;

                                        for (int intx = XX(corner1); intx <= XX(corner2); ++intx)
                                        {
                                            // Compute distance to the center of the blob
                                            // Compute blob value at that distance
                                            double d2 = A1D_ELEM(x2precalculated,intx) + y2z2;

                                            if (d2 > blobRadiusSquared)
                                                continue;
                                            int aux = (int)(d2 * iDeltaSqrt + 0.5);//Same as ROUND but avoid comparison
                                            double w = VEC_ELEM(blobTableSqrt, aux)*threadParams->weight *wModulator;

                                            // Look for the location of this logical index
                                            // in the physical layout
#ifdef DEBUG

                                            std::cout << "   gcurrent=" << gcurrent.transpose()
                                            << " d=" << d << std::endl;
                                            std::cout << "   1: intx=" << intx
                                            << " inty=" << inty
                                            << " intz=" << intz << std::endl;
#endif

                                            int ix=A1D_ELEM(xWrapped,intx);
#ifdef DEBUG

                                            std::cout << "   2: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << std::endl;
#endif

                                            bool conjugate=false;
                                            int izp, iyp, ixp;
                                            if (ix > xsize_1)
                                            {
                                                izp = izneg;
                                                iyp = iyneg;
                                                ixp = A1D_ELEM(xNegWrapped,intx);
                                                conjugate=true;
                                                fixSize = size1;
                                            }
                                            else
                                            {
                                                izp=iz;
                                                iyp=iy;
                                                ixp=ix;
                                                fixSize = size2;
                                            }
#ifdef DEBUG
                                            std::cout << "   3: ix=" << ix << " iy=" << iy
                                            << " iz=" << iz << " conj="
                                            << conjugate << std::endl;
#endif

                                            // Add the weighted coefficient
                                            if (reprocessFlag)
                                            {
                                                // Use VoutFourier as temporary to save the memory
                                                double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, izp,iyp,ixp));
                                                DIRECT_A3D_ELEM(fourierWeights, izp,iyp,ixp) += (w * ptrOut[0]);
                                            }
                                            else
                                            {
//                                              double wEffective=w*wCTF;
                                                size_t memIdx=fixSize + ixp;//YXSIZE(VoutFourier)*(izp)+((iyp)*XSIZE(VoutFourier))+(ixp);

                                                double *ptrOut=(double *)&(DIRECT_A1D_ELEM(VoutFourier, memIdx));
                                                double *ptrInVol=(double *)&(DIRECT_A1D_ELEM(parent->fftVin, memIdx));

                                                  ptrOut[0] += w * ((wCTF*ptrIn[0]-ptrInVol[0])*(wCTF*ptrIn[0]-ptrInVol[0]));
//                                                ptrOut[0] += wEffective * ((ptrIn[0]-ptrInVol[0])*(ptrIn[0]-ptrInVol[0])+(ptrIn[1]-ptrInVol[1])*(ptrIn[1]-ptrInVol[1]));
//                                                ptrOut[0] += wEffective * (ptrIn[0]);
//                                                ptrOut[0] += wEffective * std::abs(ptrIn[0] - ptrInVol[0]);
//                                                ptrOut[0] += wEffective * (ptrInVol[0]);
                                                DIRECT_A1D_ELEM(fourierWeights, memIdx) += w;
                                                if (conjugate)
//                                                	 ptrOut[1]-=wEffective * (ptrIn[1]);
                                                	ptrOut[1]+= w *((-wCTF*ptrIn[1]-ptrInVol[1])*(-wCTF*ptrIn[1]-ptrInVol[1]));
//                                                	ptrOut[1]-=wEffective * (ptrInVol[1]);
                                                else
//                                            	   ptrOut[1]+=wEffective * (ptrIn[1]);
                                            	   ptrOut[1]+= w *((wCTF*ptrIn[1]-ptrInVol[1])*(wCTF*ptrIn[1]-ptrInVol[1]));
//                                            	   ptrOut[1]+=wEffective * (ptrInVol[1]);
//                                            	   ptrOut[1]+=wEffective * std::abs(ptrIn[1]- ptrInVol[1]);

//                                                std::cout << DIRECT_A3D_ELEM(VoutFourier, izp,iyp,ixp) << " " << DIRECT_A1D_ELEM(VoutFourier, memIdx) <<  izp << " " << iyp << " " << ixp << " " << std::endl;
//                                                std::cin.get();

//                                                if ( (izp == 80) && (iyp == 80) && (ixp == 80))
//                                                	std::cout << ptrIn[0] << " " << ptrInVol[0] << " " << ptrIn[1] << " " << ptrInVol[1] << " " << w << " " << wCTF << std::endl;
                                            }
                                        }
                                    }
                                }
                            }
                    }

                    pthread_mutex_lock( &(parent->workLoadMutex) );

                    for ( int w = (minAssignedRow - minSeparation) ; w < minAssignedRow ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < (int)YSIZE(*paddedFourier) ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }

                    for ( int w = maxAssignedRow+1 ; w <= (maxAssignedRow+minSeparation) ; w ++ )
                    {
                        if ( ( w >= 0 ) && ( w < (int)YSIZE(*paddedFourier) ))
                        {
                            if ( statusArray[w] > 0 )
                            {
                                statusArray[w]--;
                            }
                        }
                    }

                    pthread_mutex_unlock( &(parent->workLoadMutex) );

                }
                while (!breakCase);
                break;
            }
        default:
            break;
        }

        barrier_wait( barrier );
    }
    while ( 1 );
}

//#define DEBUG
void ProgVolVariability::processImages( int firstImageIndex, int lastImageIndex, bool saveFSC, bool reprocessFlag)
{

	MultidimArray< std::complex<double> > *paddedFourier;

    int repaint = (int)ceil((double)SF.size()/60);

    bool processed;
    int imgno = 0;
    int imgIndex = firstImageIndex;

    // This index tells when to save work for later FSC usage
    int FSCIndex = (firstImageIndex + lastImageIndex)/2;
    // Index of the image that has just been processed. Used for
    // FSC purposes
    int current_index;

    do
    {

//PRELOAD_IMAGE
        threadOpCode = PRELOAD_IMAGE;

        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( imgIndex <= lastImageIndex )
            {
                th_args[nt].imageIndex = imgIndex;
                th_args[nt].reprocessFlag = reprocessFlag;
                imgIndex++;
            }
            else
            {
                th_args[nt].imageIndex = -1;
            }
        }

        // Awaking sleeping threads
        barrier_wait( &barrier );
        // here each thread is reading a different image and compute fft
        // Threads are working now, wait for them to finish
        // processing current projection
        barrier_wait( &barrier );

        // each threads have read a different image and now
        // all the thread will work in a different part of a single image.

//PROCESS_IMAGE
        threadOpCode = PROCESS_IMAGE;

        processed = false;

        for ( int nt = 0 ; nt < numThreads ; nt ++ )
        {
            if ( th_args[nt].read == 2 )
                processed = true;
            else if ( th_args[nt].read == 1 )
            {
                processed = true;
                if (verbose && imgno++%repaint==0)
                    progress_bar(imgno);

                double weight = th_args[nt].localweight;
                paddedFourier = th_args[nt].localPaddedFourier;
                current_index = th_args[nt].imageIndex;
                Matrix2D<double> *Ainv = th_args[nt].localAInv;

                //#define DEBUG22
#ifdef DEBUG22

                {
                    static int ii=0;
                    if(ii%1==0)
                    {
                        FourierImage save22;
                        //save22()=*paddedFourier;
                        save22().alias(*paddedFourier);
                        save22.write((std::string) integerToString(ii)  + "_padded_fourier.spi");
                    }
                    ii++;
                }
#endif
                #undef DEBUG22

                // Initialized just once
                if ( statusArray == NULL )
                {
                    statusArray = (int *) malloc ( sizeof(int) * paddedFourier->ydim );
                }

                // Determine how many rows of the fourier
                // transform are of interest for us. This is because
                // the user can avoid to explore at certain resolutions
                size_t conserveRows=(size_t)ceil( 2.0*((double)paddedFourier->ydim * maxResolution));
                conserveRows=(size_t)ceil((double)conserveRows/2.0);

                // Loop over all symmetries
                for (size_t isym = 0; isym < R_repository.size(); isym++)
                {
                    rowsProcessed = 0;

                    // Compute the coordinate axes of the symmetrized projection
                    Matrix2D<double> A_SL=R_repository[isym]*(*Ainv);

                    // Fill the thread arguments for each thread
                    for ( int th = 0 ; th < numThreads ; th ++ ) // JV I think this for is incorrect all the threads have the same A_SL matrix
                    {
                        // Passing parameters to each thread
                        th_args[th].symmetry = &A_SL;
                        th_args[th].paddedFourier = paddedFourier;
                        th_args[th].weight = weight;
                        th_args[th].reprocessFlag = reprocessFlag;
                    }

                    // Init status array
                    for (size_t i = 0 ; i < paddedFourier->ydim ; i ++ )
                    {
                        if ( i >= conserveRows && i < (paddedFourier->ydim-conserveRows))
                        {
                            // -2 means "discarded"
                            statusArray[i] = -2;
                            rowsProcessed++;
                        }
                        else
                        {
                            statusArray[i] = 0;
                        }
                    }

                    // Awaking sleeping threads
                    barrier_wait( &barrier );
                    // Threads are working now, wait for them to finish
                    // processing current projection
                    barrier_wait( &barrier );
//PROCESS_IMAGE

                    //#define DEBUG2
#ifdef DEBUG2

                    {
                        static int ii=0;
                        if(ii%1==0)
                        {
                            Image<double> save;
                            save().alias( FourierWeights );
                            save.write((std::string) integerToString(ii)  + "_1_Weights.vol");

                            Image< std::complex<double> > save2;
                            save2().alias( VoutFourier );
                            save2.write((std::string) integerToString(ii)  + "_1_Fourier.vol");
                        }
                        ii++;
                    }
#endif
                    #undef DEBUG2

                }

                if ( current_index == FSCIndex && saveFSC )
                {
                    // Save Current Fourier, Reconstruction and Weights
                    Image<double> save;
                    save().alias( FourierWeights );
                    save.write((std::string)fn_fsc + "_1_Weights.vol");

                    Image< std::complex<double> > save2;
                    save2().alias( VoutFourier );
                    save2.write((std::string) fn_fsc + "_1_Fourier.vol");

                    finishComputations(FileName((std::string) fn_fsc + "_1_recons.vol"));
                    Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
                    transformerVol.setReal(Vout());
                    Vout().clear();
                    transformerVol.getFourierAlias(VoutFourier);
                    FourierWeights.initZeros(VoutFourier);
                    VoutFourier.initZeros();
                }
            }
        }
    }
    while ( processed );

    if( saveFSC )
    {
        // Save Current Fourier, Reconstruction and Weights
        Image<double> auxVolume;
        auxVolume().alias( FourierWeights );
        auxVolume.write((std::string)fn_fsc + "_2_Weights.vol");

        Image< std::complex<double> > auxFourierVolume;
        auxFourierVolume().alias( VoutFourier );
        auxFourierVolume.write((std::string) fn_fsc + "_2_Fourier.vol");

        finishComputations(FileName((std::string) fn_fsc + "_2_recons.vol"));

        Vout().initZeros(volPadSizeZ, volPadSizeY, volPadSizeX);
        transformerVol.setReal(Vout());
        Vout().clear();
        transformerVol.getFourierAlias(VoutFourier);
        FourierWeights.initZeros(VoutFourier);
        VoutFourier.initZeros();

        auxVolume.sumWithFile(fn_fsc + "_1_Weights.vol");
        auxVolume.sumWithFile(fn_fsc + "_2_Weights.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_1_Fourier.vol");
        auxFourierVolume.sumWithFile(fn_fsc + "_2_Fourier.vol");
        remove((fn_fsc + "_1_Weights.vol").c_str());
        remove((fn_fsc + "_2_Weights.vol").c_str());
        remove((fn_fsc + "_1_Fourier.vol").c_str());
        remove((fn_fsc + "_2_Fourier.vol").c_str());

        /*
        //Save SUM
                                    //this is an image but not an xmipp image
                                    auxFourierVolume.write((std::string)fn_fsc + "_all_Fourier.vol",
                                            false,VDOUBLE);
                                    auxVolume.write((std::string)fn_fsc + "_all_Weight.vol",
                                            false,VDOUBLE);
        //
        */
    }
}

void ProgVolVariability::correctWeight()
{
    // If NiterWeight=0 then set the weights to one
	forceWeightSymmetry(FourierWeights);
    if (NiterWeight==0)
    {
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FourierWeights)
        DIRECT_A3D_ELEM(FourierWeights, k,i,j)=1;
    }
    else
    {
        // Temporary save the Fourier of the volume
        MultidimArray< std::complex<double> > VoutFourierTmp;
        VoutFourierTmp=VoutFourier;
        forceWeightSymmetry(FourierWeights);
        // Prepare the VoutFourier
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
        {
            double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));

            if (fabs(A3D_ELEM(FourierWeights,k,i,j))>1e-3)
                ptrOut[0] = 1.0/DIRECT_A3D_ELEM(FourierWeights, k,i,j);

        }

        for (int i=1;i<NiterWeight;i++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(FourierWeights)
            A3D_ELEM(FourierWeights,k,i,j)=0;
            processImages(0, SF.size() - 1, !fn_fsc.empty(), true);
            forceWeightSymmetry(FourierWeights);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
            {
                double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
                if (fabs(A3D_ELEM(FourierWeights,k,i,j))>1e-3)
                    ptrOut[0] /= A3D_ELEM(FourierWeights,k,i,j);

            }
        }
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(VoutFourier)
        {
            // Put back the weights to FourierWeights from temporary variable VoutFourier
            double *ptrOut=(double *)&(DIRECT_A3D_ELEM(VoutFourier, k,i,j));
            A3D_ELEM(FourierWeights,k,i,j) = ptrOut[0];
        }
        VoutFourier = VoutFourierTmp;
    }
}

void ProgVolVariability::finishComputations( const FileName &out_name )
{

#ifdef DEBUG_VOL
    {
    	Image<double> save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "Weights.vol");

        Image< std::complex<double> > save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "FourierVol.vol");
    }
#endif

    // Enforce symmetry in the Fourier values as well as the weights
    // Sjors 19aug10 enforceHermitianSymmetry first checks ndim...


    //~JV
    //Read the input average map
    Image<double> Vin;
    Vin.read(fn_vol);
    Vin().setXmippOrigin();
    double corr2D_3D=pow(padding_factor_proj,2.)/
                     (imgSize* pow(padding_factor_vol,3.));


    Vout().initZeros(ZSIZE(Vin()),YSIZE(Vin()),XSIZE(Vin()));
    Vout().setXmippOrigin();

    MultidimArray< double > Vmc;
    Vmc.initZeros(volPadSizeZ,volPadSizeY,volPadSizeX);

    MultidimArray< double > Vintemp;
    Vintemp.initZeros(ZSIZE(Vin()),YSIZE(Vin()),XSIZE(Vin()));
    Vintemp.setXmippOrigin();

    transformerVol.setReal(Vmc);
    transformerVol.enforceHermitianSymmetry();

	if (fn_mask != "")
	{
		mask.read(fn_mask);
		mask().setXmippOrigin();
	}

    //forceWeightSymmetry(preFourierWeights);

    // Tell threads what to do
    //#define DEBUG_VOL1
#ifdef DEBUG_VOL1

    {
        Image<double> save;
        save().alias( FourierWeights );
        save.write((std::string) fn_out + "hermiticWeights.vol");

        Image< std::complex<double> > save2;
        save2().alias( VoutFourier );
        save2.write((std::string) fn_out + "hermiticFourierVol.vol");
    }
#endif

    //PROCESS_WEIGHTS
    threadOpCode = PROCESS_WEIGHTS;
    // Awake threads
    barrier_wait( &barrier );
    // Threads are working now, wait for them to finish
    barrier_wait( &barrier );

//~JV

    //MONTE CARLO SIMULATION
    MultidimArray< std::complex<double> > VoutFourierTmp2;
    VoutFourierTmp2 = VoutFourier;
    std::complex<double>  mean = 0;
    std::complex<double> stddev = 0;
    std::complex<double> result = 0;

    std::cout << std::endl;
    double error = 0.0;
    double mean_err = 1000;

    int it = 1;
    std::cout << "Monte Carlo simulation: " << std::endl;

    while( it <= NiterMC )
    {

		FOR_ALL_ELEMENTS_IN_ARRAY3D(VoutFourierTmp2)  // Simulated Volume
		{
			mean = A3D_ELEM(fftVin, k, i, j)*corr2D_3D;
			stddev = (A3D_ELEM(VoutFourierTmp2, k, i, j)); // We compute the std from the variance

			stddev.real( std::sqrt( std::abs(stddev.real()))*corr2D_3D );
			stddev.imag( std::sqrt( std::abs(stddev.imag()))*corr2D_3D );

			rand_normal(((double*) &mean)[0], ((double*) &stddev)[0],
					((double*) &result)[0]);
			rand_normal(((double*) &mean)[1], ((double*) &stddev)[1],
					((double*) &result)[1]);
			A3D_ELEM(VoutFourier,k,i,j) = result;//result; // A3D_ELEM(fftVin,k,i,j);
		} 		// Simulated Volume


//Fourier Transform of simulated volume
	    transformerVol.inverseFourierTransform();
	    CenterFFT(Vmc,false);

	    // Correct by the Fourier transform of the blob
	    Vmc.setXmippOrigin();

	    //Remove the padding
	    Vmc.selfWindow(FIRST_XMIPP_INDEX(imgSize),FIRST_XMIPP_INDEX(imgSize),
	                      FIRST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize),
	                      LAST_XMIPP_INDEX(imgSize),LAST_XMIPP_INDEX(imgSize));

	    double pad_relation= ((double)padding_factor_proj/padding_factor_vol);
	    pad_relation = (pad_relation * pad_relation * pad_relation);

	    MultidimArray<double> &mVout=Vmc;
	    double ipad_relation=1.0/pad_relation;
	    double meanFactor2=0;
	    FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
	    {

	        double radius=sqrt((double)(k*k+i*i+j*j));
	        double aux=radius*iDeltaFourier;
	        double factor = Fourier_blob_table(ROUND(aux));
	        double factor2=(pow(Sinc(radius/(2*(imgSize))),2));
	        if (NiterWeight!=0)
	        {
	            A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor2*factor);
	            meanFactor2+=factor2;
	        }
	        else
	            A3D_ELEM(mVout,k,i,j) /= (ipad_relation*factor);
	    }

	    if (NiterWeight!=0) //Esto esta mal tiene que estar dentro de un bucle xq aqui k,i,j no tienen sentido
	    {
	        meanFactor2/=MULTIDIM_SIZE(mVout);
	        FOR_ALL_ELEMENTS_IN_ARRAY3D(mVout)
	        	A3D_ELEM(mVout,k,i,j) *= meanFactor2;
	    }

	    int numElem = 0;
	    FOR_ALL_ELEMENTS_IN_ARRAY3D(Vmc) //Esto deberia juntarlo con el bucle superior
	    {
	    	//This is variance
	    	A3D_ELEM(Vout(),k,i,j) += ((A3D_ELEM(Vmc,k,i,j)-A3D_ELEM(Vin(),k,i,j))*(A3D_ELEM(Vmc,k,i,j)-A3D_ELEM(Vin(),k,i,j)));
	    	A3D_ELEM(Vintemp,k,i,j) += (A3D_ELEM(Vmc,k,i,j));

			if (fn_mask != "") //incluir en el if las lineas anteriores para acelear el procesamiento
	    		if ( A3D_ELEM(mask(),k,i,j) != 0 )
	    		{
	    			error += (std::abs(A3D_ELEM(Vintemp,k,i,j)/(double)it - A3D_ELEM(Vin(),k,i,j)));
	    	    	numElem++;
	    		}
	    		else
	    			;
	    	else
	    	{
    			error += (std::abs(A3D_ELEM(Vintemp,k,i,j)/(double)it - A3D_ELEM(Vin(),k,i,j)));
    	    	numElem++;
	    	}
	    }

	    mean_err = (error/(numElem))*100.0;
	    std::cout << "Iteration MC : " << it << " Max iter: " << NiterMC << " Error : " << mean_err << std::endl;
	    Vmc.initZeros(volPadSizeZ,volPadSizeY,volPadSizeX); //we want to reuse the Vin() memory
	    Vmc.setXmippOrigin();
	    error = 0.0;
	    it=it+1;
    }
    //~JV

    Vout()  /= (double)it;
    Vintemp /= (double)it;

    FileName fn_variance = fn_out.removeFileFormat().removeLastExtension()+"_variance.vol";
    Vout.write(fn_variance);

    fn_variance = fn_out.removeFileFormat().removeLastExtension()+ "_mean.vol";
    Vin() = Vintemp;
    Vin.write(fn_variance);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(Vout())
    	A3D_ELEM(Vout(),k,i,j) = std::sqrt(A3D_ELEM(Vout(),k,i,j));

    fn_variance = fn_out.removeFileFormat().removeLastExtension()+ "_std.vol";
    Vout.write(fn_variance);


    FOR_ALL_ELEMENTS_IN_ARRAY3D(Vintemp)
    	A3D_ELEM(Vout(),k,i,j) = std::abs(A3D_ELEM(Vintemp,k,i,j))/A3D_ELEM(Vout(),k,i,j);

    fn_variance = fn_out.removeFileFormat().removeLastExtension()+ "_snr.vol";
    Vout.write(fn_variance);

    std::cout << std::endl;


    #ifdef DEBUG_VOL2 //poner estas lineas
    {
    	Vintemp /=(double)it;
        Image<double> save;
        save().alias( Vintemp );
        save.write((std::string) fn_out + "Vintemp.vol");
        Vin.write((std::string) fn_out + "Vin.vol");
    }
    #endif

}

void ProgVolVariability::setIO(const FileName &fn_in, const FileName &fn_out)
{
    this->fn_sel = fn_in;
    this->fn_out = fn_out;
}

void ProgVolVariability::forceWeightSymmetry(MultidimArray<double> &FourierWeights)
{
    int yHalf=YSIZE(FourierWeights)/2;
    if (YSIZE(FourierWeights)%2==0)
        yHalf--;
    int zHalf=ZSIZE(FourierWeights)/2;
    if (ZSIZE(FourierWeights)%2==0)
        zHalf--;
    int zsize=(int)ZSIZE(FourierWeights);
    int zsize_1=zsize-1;
    int ysize_1=(int)YSIZE(FourierWeights)-1;
    for (int k=0; k<zsize; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        for (int i=1; i<=yHalf; i++)
        {
            int isym=intWRAP(-i,0,ysize_1);
            double mean=0.5*(
                            DIRECT_A3D_ELEM(FourierWeights,k,i,0)+
                            DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0));
            DIRECT_A3D_ELEM(FourierWeights,k,i,0)=
                DIRECT_A3D_ELEM(FourierWeights,ksym,isym,0)=mean;
        }
    }
    for (int k=1; k<=zHalf; k++)
    {
        int ksym=intWRAP(-k,0,zsize_1);
        double mean=0.5*(
                        DIRECT_A3D_ELEM(FourierWeights,k,0,0)+
                        DIRECT_A3D_ELEM(FourierWeights,ksym,0,0));
        DIRECT_A3D_ELEM(FourierWeights,k,0,0)=
            DIRECT_A3D_ELEM(FourierWeights,ksym,0,0)=mean;
    }
}


void ProgVolVariability::rand_normal(double& mean, double& stddev, double& result)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;

    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            result = n1*stddev + mean;
            n2_cached = 1;
        }
    }
    else
    {
        n2_cached = 0;
        result = n2*stddev + mean;
    }
}


