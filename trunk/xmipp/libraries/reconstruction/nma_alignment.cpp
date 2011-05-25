/***************************************************************************
 *
 * Authors:    Slavica Jonic                slavica.jonic@impmc.jussieu.fr
 *             Carlos Oscar Sanchez Sorzano coss.eps@ceu.es
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
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#include "nma_alignment.h"
#include "data/metadata_extension.h"
#include "program_extension.h"

#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include "../../external/condor/Solver.h"
#include "../../external/condor/tools.h"

// Empty constructor =======================================================
ProgNmaAlignment::ProgNmaAlignment()
{
    rangen = 0;
    resume=false;
    currentImgName="";
    each_image_produces_an_output = false;
    produces_an_output = true;
    progVolumeFromPDB = new ProgPdbConverter();
}

ProgNmaAlignment::~ProgNmaAlignment()
{
	delete progVolumeFromPDB;
}

// Params definition ============================================================
void ProgNmaAlignment::defineParams()
{
	addUsageLine("Compute deformation parameters according to a set of NMA modes");
	defaultComments["-o"].clear();
	defaultComments["-o"].addComment("Metadata with output alignment and deformations");
    XmippMetadataProgram::defineParams();
    addParamsLine("   --pdb <PDB_filename>                : PDB Model to compute NMA");
    addParamsLine("  [--resume]                           : Resume processing");
    addParamsLine("==Generation of the deformed volumes==");
    addParamsLine("   --modes <filename>                  : File with a list of mode filenames");
    addParamsLine("  [--deformation_scale <s=1>]          : Scaling factor to scale NMA deformation amplitudes");
    addParamsLine("  [--sampling_rate <Ts=1>]             : in Angstroms/pixel");
    addParamsLine("  [--filterVol <cutoff=15.>]           : Filter the volume after deforming. Default cut-off is 15 A.");
    addParamsLine("  [--centerPDB]                        : Center the PDB structure");
    addParamsLine("  [--fixed_Gaussian <std=-1>]          : For pseudo atoms fixed_Gaussian must be used.");
    addParamsLine("                                       : Default standard deviation <std> is read from PDB file.");
    addParamsLine("==Angular assignment and mode detection==");
    addParamsLine("  [--mask <m=\"\">]                    : 2D Mask applied to the reference images of the deformed volume");
    addParamsLine("  [--gaussian_Fourier <s=0.5>]         : Weighting sigma in Fourier space");
    addParamsLine("  [--gaussian_Real    <s=0.5>]         : Weighting sigma in Real space");
    addParamsLine("  [--zerofreq_weight  <s=0.>]          : Zero-frequency weight");
    addExampleLine("xmipp_nma_alignment -i images.sel --pdb 2tbv.pdb --modes modelist.txt --deformation_scale 1000 --sampling_rate 6.4 -o output.txt --resume");
}

// Read arguments ==========================================================
void ProgNmaAlignment::readParams()
{
    XmippMetadataProgram::readParams();
    fnPDB = getParam("--pdb");
    fnModeList = getParam("--modes");
    resume = checkParam("--resume");
    scale_defamp = getDoubleParam("--deformation_scale");
    sampling_rate = getDoubleParam("--sampling_rate");
    fnmask = getParam("--mask");
    gaussian_DFT_sigma = getDoubleParam( "--gaussian_Fourier");
    gaussian_Real_sigma = getDoubleParam( "--gaussian_Real");
    weight_zero_freq = getDoubleParam( "--zerofreq_weight");
    do_centerPDB = checkParam("--centerPDB");
    do_FilterPDBVol = checkParam("--filterVol");
    if (do_FilterPDBVol)
        cutoff_LPfilter = getDoubleParam("--filterVol");
    useFixedGaussian = checkParam("--fixed_Gaussian");
    if (useFixedGaussian)
        sigmaGaussian = getDoubleParam("--fixed_Gaussian");

}

// Show ====================================================================
void ProgNmaAlignment::show()
{
    XmippMetadataProgram::show();
    std::cout
    << "PDB:                 " << fnPDB               << std::endl
    << "Resume:              " << resume              << std::endl
    << "Mode list:           " << fnModeList          << std::endl
    << "Amplitude scale:     " << scale_defamp        << std::endl
    << "Sampling rate:       " << sampling_rate       << std::endl
    << "Mask:                " << fnmask              << std::endl
    << "Gaussian Fourier:    " << gaussian_DFT_sigma  << std::endl
    << "Gaussian Real:       " << gaussian_Real_sigma << std::endl
    << "Zero-frequency weight:"<< weight_zero_freq    << std::endl
    << "Center PDB:          " << do_centerPDB        << std::endl
    << "Filter PDB volume    " << do_FilterPDBVol     << std::endl
    << "Use fixed Gaussian:  " << useFixedGaussian    << std::endl
    << "Sigma of Gaussian:   " << sigmaGaussian       << std::endl
    ;
}

// Produce side information ================================================
ProgNmaAlignment *global_NMA_prog;

void ProgNmaAlignment::createWorkFiles()
{
    MetaData mdTodo, mdDone;
    mdTodo = mdIn;
    if (exists("nmaDone.xmd") && resume)
    {
        mdDone.read("nmaDone.xmd");
        mdTodo.subtraction(mdDone, MDL_IMAGE);
    }
    else//if not exists create metadata only with headers
    {
        mdDone.addLabel(MDL_IMAGE);
        mdDone.addLabel(MDL_ENABLED);
        mdDone.addLabel(MDL_ANGLEROT);
        mdDone.addLabel(MDL_ANGLETILT);
        mdDone.addLabel(MDL_ANGLEPSI);
        mdDone.addLabel(MDL_SHIFTX);
        mdDone.addLabel(MDL_SHIFTY);
        mdDone.addLabel(MDL_NMA);
        mdDone.addLabel(MDL_COST);
        mdDone.write("nmaDone.xmd");
    }
    mdIn = mdTodo;
}

void ProgNmaAlignment::preProcess()
{
	MetaData SF(fnModeList);
	numberOfModes=SF.size();
    // Get the size of the images in the selfile
    int ydim, zdim;
    size_t ndim;
    ImgSize(mdIn, imgSize, ydim, zdim, ndim);
    // Set the pointer of the program to this object
    global_NMA_prog = this;
    //create some neededs files
    createWorkFiles();
}

void ProgNmaAlignment::finishProcessing()
{
	XmippMetadataProgram::finishProcessing();
	rename("nmaDone.xmd",fn_out.c_str());
}

void runSystem(const String &program, const String &arguments, bool useSystem = true)
{
    if (useSystem)
    {
        String cmd = formatString("%s %s", program.c_str(), arguments.c_str());
#ifdef DEBUG
        std::cerr << std::endl << ">>> RUNNING EXTERNALLY: " << cmd << std::endl;
#endif
        system(cmd.c_str());
    }
    else
    {
#ifdef DEBUG
    	std::cerr << std::endl << ">>> RUNNING INTERNALLY: " << program << std::endl;
#endif
        runProgram(program, arguments);
    }
}

// Create deformed PDB =====================================================
FileName ProgNmaAlignment::createDeformedPDB(int pyramidLevel) const
{
    String program;
    String arguments;
    FileName fnRandom;
    fnRandom.initUniqueName(nameTemplate);
    const char * randStr = fnRandom.c_str();

    program = "xmipp_pdb_nma_deform";
    arguments = formatString("--pdb %s -o deformedPDB_%s.pdb --nma %s --deformations ",
                             fnPDB.c_str(), randStr, fnModeList.c_str());
    for (int i=0; i<VEC_XSIZE(trial)-5; ++i)
        arguments += floatToString(trial(i)*scale_defamp)+" ";
    runSystem(program, arguments, false);

    program = "xmipp_volume_from_pdb";
    arguments = formatString("-i deformedPDB_%s.pdb --size %i --sampling %f -v 0",
                             randStr, imgSize, sampling_rate);

    if (do_centerPDB)
        arguments.append(" --centerPDB ");

    if (useFixedGaussian)
    {
        arguments.append(" --fixed_Gaussian ");
        if (sigmaGaussian >= 0)
            arguments += formatString("%f --intensityColumn Bfactor", sigmaGaussian);
    }
    progVolumeFromPDB->read(arguments);
    progVolumeFromPDB->tryRun();

    if (do_FilterPDBVol)
    {
        program = "xmipp_transform_filter";
        arguments = formatString("-i deformedPDB_%s.vol --sampling %f --fourier low_pass %f raised_cosine 0.1 -v 0",
                                 randStr, sampling_rate, cutoff_LPfilter);
        runSystem(program, arguments, false);
    }

    if (pyramidLevel != 0)
    {
        Image<double> I;
        FileName fnDeformed = formatString("deformedPDB_%s.vol", randStr);
        I.read(fnDeformed);
        selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
        I.write(fnDeformed);
    }

    return fnRandom;
}

// Perform complete search =================================================
void ProgNmaAlignment::performCompleteSearch(
    const FileName &fnRandom, int pyramidLevel) const
{
    String program;
    String arguments;
    const char * randStr = fnRandom.c_str();

    // Reduce the image
    FileName fnDown = formatString("downimg_%s.xmp", fnRandom.c_str());
    if (pyramidLevel != 0)
    {
        Image<double> I;
        I.read(currentImgName);
        selfPyramidReduce(BSPLINE3, I(), pyramidLevel);
        I.write(fnDown);
    }
    else
        link(currentImgName.c_str(),fnDown.c_str());

    mkdir(((std::string)"ref" + fnRandom).c_str(),S_IRWXU);

    program = "xmipp_angular_project_library";
    arguments = formatString("-i deformedPDB_%s.vol -o ref%s/ref%s.stk --sampling_rate 25 -v 0",
                             randStr, randStr, randStr);
    runSystem(program, arguments, false);

    const char * refSelStr = formatString("ref%s/ref%s.doc", randStr, randStr).c_str();

    if (fnmask != "")
    {
        program = "xmipp_transform_mask";
        arguments = formatString("-i %s --mask binary_file %s", refSelStr, fnmask.c_str());
        runSystem(program, arguments, false);
    }

    // Perform alignment
    program = "xmipp_angular_discrete_assign";
    arguments = formatString("-i downimg_%s.xmp --ref %s -o angledisc_%s.txt --psi_step 5 --max_shift_change %i --search5D -v 0",
        randStr, refSelStr, randStr, ROUND((double)imgSize/(10.0*pow(2.0,(double)pyramidLevel))));
    runSystem(program, arguments,false);
}

// Continuous assignment ===================================================
double ProgNmaAlignment::performContinuousAssignment(
    const FileName &fnRandom, int pyramidLevel) const
{
    // Perform alignment
    const char * randStr = fnRandom.c_str();
    String program = "xmipp_angular_continuous_assign";
    String arguments = formatString("-i angledisc_%s.txt --ref deformedPDB_%s.vol -o anglecont_%s.txt --gaussian_Fourier %f --gaussian_Real %f --zerofreq_weight %f -v 0",
                                    randStr, randStr, randStr, gaussian_DFT_sigma, gaussian_Real_sigma, weight_zero_freq);
    runSystem(program, arguments, false);

    // Pick up results
    MetaData DF;
    DF.read("anglecont_"+fnRandom+".txt");
    size_t objId = DF.firstObject();
    DF.getValue(MDL_ANGLEROT,trial(VEC_XSIZE(trial)-5),objId);
    DF.getValue(MDL_ANGLETILT,trial(VEC_XSIZE(trial)-4),objId);
    DF.getValue(MDL_ANGLEPSI,trial(VEC_XSIZE(trial)-3),objId);
    DF.getValue(MDL_SHIFTX,trial(VEC_XSIZE(trial)-2),objId);
    trial(VEC_XSIZE(trial)-2)*=pow(2.0,(double)pyramidLevel);
    DF.getValue(MDL_SHIFTY,trial(VEC_XSIZE(trial)-1),objId);
    trial(VEC_XSIZE(trial)-1)*=pow(2.0,(double)pyramidLevel);
    double tempvar;
    DF.getValue(MDL_COST,tempvar,objId);
    return tempvar;
}

void ProgNmaAlignment::updateBestFit(double fitness,int dim)
{
    if (fitness < fitness_min(0))
    {
        fitness_min(0) = fitness;
        trial_best = trial;
    }
}

// Compute fitness =========================================================
double ObjFunc_nma_alignment::eval(Vector X, int *nerror)
{
    int dim = global_NMA_prog->numberOfModes;

    for (int i=0; i<dim; i++)
    {
        global_NMA_prog->trial(i) = X[i];
    }

    int pyramidLevelDisc = 1;
    int pyramidLevelCont = (global_NMA_prog->currentStage == 1) ? 1 : 0;

    FileName fnRandom = global_NMA_prog->createDeformedPDB(pyramidLevelCont);
    const char * randStr = fnRandom.c_str();

    if (global_NMA_prog->currentStage==1)
    {
        global_NMA_prog->performCompleteSearch(fnRandom,pyramidLevelDisc);
    }
    else
    {
        double rot, tilt, psi, xshift, yshift;
        MetaData DF;

        rot = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-5);
        tilt = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-4);
        psi = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-3);
        xshift = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-2);
        yshift = global_NMA_prog->bestStage1(VEC_XSIZE(global_NMA_prog->bestStage1)-1);

        size_t objId=DF.addObject();
        FileName fnDown = formatString("downimg_%s.xmp", randStr);
        DF.setValue(MDL_IMAGE,fnDown,objId);
        DF.setValue(MDL_ENABLED,1,objId);
        DF.setValue(MDL_ANGLEROT,rot,objId);
        DF.setValue(MDL_ANGLETILT,tilt,objId);
        DF.setValue(MDL_ANGLEPSI,psi,objId);
        DF.setValue(MDL_SHIFTX,xshift,objId);
        DF.setValue(MDL_SHIFTY,yshift,objId);

        DF.write(formatString("angledisc_%s.txt", randStr));
        link(global_NMA_prog->currentImgName.c_str(), fnDown.c_str());
    }
    double fitness = global_NMA_prog->performContinuousAssignment(fnRandom, pyramidLevelCont);

    runSystem("rm", formatString("-rf *%s* &", randStr));

    global_NMA_prog->updateBestFit(fitness,dim);

    return fitness;
}

ObjFunc_nma_alignment::ObjFunc_nma_alignment(int _t, int _n)
{}

void ProgNmaAlignment::processImage(const FileName &fnImg, const FileName &fnImgOut, size_t objId)
{
    double rhoStart=1e-0, rhoEnd=1e-3;

    int niter=1000;

    ObjectiveFunction *of;

    int dim=numberOfModes;

    parameters.initZeros(dim+5);
    currentImgName = fnImg;
    sprintf(nameTemplate, "_node%d_img%ld_XXXXXX", rangen, objId);

    trial.initZeros(dim+5);
    trial_best.initZeros(dim+5);

    fitness_min.initZeros(1);
    fitness_min(0) = 1000000.0;

    currentStage=1;
#ifdef DEBUG
    std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
    <<" processing image " << fnImg <<"(" << objId << ")"
    << " at stage: " << currentStage << std::endl;
#endif
    of=new ObjFunc_nma_alignment(1,dim);

    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++)
        of->xStart[i]=0.;

#ifdef DEBUG
    strcpy(of->name,("OF1_"+integerToString(rangen)).c_str());
    of->setSaveFile();
#endif

    CONDOR(rhoStart, rhoEnd, niter, of);
#ifdef DEBUG
    of->printStats();
    FILE *ff = fopen(("res1_"+integerToString(rangen)+".txt").c_str(),"w");
    fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
    fclose(ff);
#endif

    double fitness=of->valueBest;
    double *dd=of->xBest;

    bestStage1 = trial = parameters = trial_best;

    delete of;

    currentStage = 2;
#ifdef DEBUG
    std::cerr << std::endl << "DEBUG: ===== Node: " << rangen
    <<" processing image " << fnImg <<"(" << objId << ")"
    << " at stage: " << currentStage << std::endl;
#endif

    fitness_min(0) = 1000000.0;

    of=new ObjFunc_nma_alignment(1,dim);

    of->xStart.setSize(dim);
    for (int i=0; i<dim; i++)
        of->xStart[i]=parameters(i);
#ifdef DEBUG
    strcpy(of->name,("OF2_"+integerToString(rangen)).c_str());
    of->setSaveFile();
#endif

    rhoStart=1e-3, rhoEnd=1e-4;
    CONDOR(rhoStart, rhoEnd, niter, of);
#ifdef DEBUG
    of->printStats();
    ff=fopen(("res2_"+integerToString(rangen)+".txt").c_str(),"w");
    fprintf(ff,"%s & %i & %i & (%i) & %e \\\\\n", of->name, of->dim(), of->getNFE(), of->getNFE2(), of->valueBest);
    fclose(ff);
#endif

    fitness=of->valueBest;
    dd=of->xBest;
#ifdef DEBUG
    std::cout << "Best fitness = " << fitness << std::endl;
    for (int i=0; i<dim; i++)
    {
        std::cout << "Best deformations = " << dd[i] << std::endl;
    }
#endif

    trial = trial_best;

    for (int i=dim; i<dim+5; i++)
    {
        parameters(i-dim)=trial_best(i);
    }

    for (int i=0; i<dim; i++)
    {
        parameters(5+i)=trial_best(i)*scale_defamp;
    }

    parameters.resize(VEC_XSIZE(parameters)+1);
    parameters(VEC_XSIZE(parameters)-1)=fitness_min(0);

    writeImageParameters(fnImg);
    delete of;
}

void ProgNmaAlignment::writeImageParameters(const FileName &fnImg)
{
    MetaData md;
    size_t objId = md.addObject();
    md.setValue(MDL_IMAGE,fnImg,objId);
    md.setValue(MDL_ENABLED,1,objId);
    md.setValue(MDL_ANGLEROT,parameters(0),objId);
    md.setValue(MDL_ANGLETILT,parameters(1),objId);
    md.setValue(MDL_ANGLEPSI,parameters(2),objId);
    md.setValue(MDL_SHIFTX,parameters(3),objId);
    md.setValue(MDL_SHIFTY,parameters(4),objId);

    int dim = numberOfModes;
    std::vector<double> vectortemp;

    for (int j = 5; j < 5+dim; j++)
    {
        vectortemp.push_back(parameters(j));
    }

    md.setValue(MDL_NMA,vectortemp,objId);
    md.setValue(MDL_COST,parameters(5+dim),objId);

    md.append("nmaDone.xmd");
}
