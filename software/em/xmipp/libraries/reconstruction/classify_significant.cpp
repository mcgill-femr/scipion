/***************************************************************************
 *
 * Authors:    Carlos Oscar Sanchez Sorzano coss@cnb.csic.es
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

#include "classify_significant.h"
#include <data/mask.h>

// Empty constructor =======================================================
ProgClassifySignificant::~ProgClassifySignificant()
{
	for (size_t i=0; i<projector.size(); i++)
		delete projector[i];
	for (size_t i=0; i<subsetProjections.size(); i++)
		delete subsetProjections[i];
	for (size_t i=0; i<Iexp.size(); i++)
		delete Iexp[i];
}

// Read arguments ==========================================================
void ProgClassifySignificant::readParams()
{
    fnVols = getParam("--ref");
    fnIds = getParam("--id");
    fnAngles = getParam("--angles");
    fnOut = getParam("-o");
    pad = getIntParam("--padding");
    wmin = getDoubleParam("--minWeight");
    onlyIntersection = checkParam("--onlyIntersection");
    numVotes = getIntParam("--votes");
}

// Show ====================================================================
void ProgClassifySignificant::show()
{
    if (!verbose)
        return;
    std::cout
    << "Reference volumes:   " << fnVols            << std::endl
    << "Ids:                 " << fnIds             << std::endl
    << "Angles:              " << fnAngles          << std::endl
    << "Output:              " << fnOut             << std::endl
    << "Padding factor:      " << pad               << std::endl
	<< "Min. Weight:         " << wmin              << std::endl
    ;
}

// usage ===================================================================
void ProgClassifySignificant::defineParams()
{
    addUsageLine("Classify a set of images into different classes. See protocol_reconstruct_heterogeneous");
    addParamsLine("   --ref <metadata>            : Reference volumes");
    addParamsLine("   --id <metadata>             : List of itemIds to classified. Sorted.");
    addParamsLine("   --angles <metadata>         : Angles assigned. Each image should have one or several angles");
    addParamsLine("                               : for each volume. The assignment per volume should be organized");
    addParamsLine("                               : in different blocks");
    addParamsLine("   -o <metadata>               : Output metadata with a set of angles per volume");
    addParamsLine("  [--votes <numVotes=5>]       : Minimum number of votes to consider an image belonging to a volume");
    addParamsLine("  [--onlyIntersection]         : Flag to select only the images belonging only to the set intersection");
    addParamsLine("  [--padding <p=2>]            : Padding factor");
    addParamsLine("  [--minWeight <w=0.1>]        : Minimum weight");
}

// Produce side information ================================================
void ProgClassifySignificant::produceSideInfo()
{
	if (verbose>0)
		std::cerr << "Producing side info ..." << std::endl;
    // Read the reference volumes
    Image<double> V;
    MetaData mdVols, mdAngles, mdAnglesSorted;
    mdVols.read(fnVols);
    FileName fnVol;
    int i=1;
    FOR_ALL_OBJECTS_IN_METADATA(mdVols)
    {
    	mdVols.getValue(MDL_IMAGE,fnVol,__iter.objId);
    	std::cout << fnVol << std::endl;
        V.read(fnVol);
        V().setXmippOrigin();
        projector.push_back(new FourierProjector(V(),pad,0.5,BSPLINE3));
        currentRowIdx.push_back(0);

        mdAngles.read(formatString("angles_%02d@%s",i,fnAngles.c_str()));
        mdAnglesSorted.sort(mdAngles, MDL_ITEM_ID, true);
        VMetaData *vmd=new VMetaData();
        mdAnglesSorted.asVMetaData(*vmd);
        setAngles.push_back(*vmd);
        classifiedAngles.push_back(*(new VMetaData()));

        subsetAngles.push_back(*(new VMetaData()));
        subsetProjectionIdx.push_back(* (new std::vector<size_t>));
        i+=1;
    }

    // Read the Ids
    MetaData mdIds;
    mdIds.read(fnIds);
    mdIds.getColumnValues(MDL_PARTICLE_ID,setIds);
}

//#define DEBUG
void ProgClassifySignificant::generateProjection(size_t volumeIdx, size_t poolIdx, MDRow &currentRow)
{
	double rot, tilt, psi, x, y;
	bool flip;
	Matrix2D<double> A;
	Image<double> &Iaux = *Iexp[0];
	int xdim = (int)XSIZE(Iaux());

	currentRow.getValue(MDL_ANGLE_ROT,rot);
	currentRow.getValue(MDL_ANGLE_TILT,tilt);
	currentRow.getValue(MDL_ANGLE_PSI,psi);
	currentRow.getValue(MDL_SHIFT_X,x);
	currentRow.getValue(MDL_SHIFT_Y,y);
	currentRow.getValue(MDL_FLIP,flip);
	A.initIdentity(3);
	MAT_ELEM(A,0,2)=x;
	MAT_ELEM(A,1,2)=y;
	if (flip)
	{
		MAT_ELEM(A,0,0)*=-1;
		MAT_ELEM(A,0,1)*=-1;
		MAT_ELEM(A,0,2)*=-1;
	}
	projectVolume(*(projector[volumeIdx]), Paux, xdim, xdim,  rot, tilt, psi);

	if (poolIdx>=subsetProjections.size())
		subsetProjections.push_back(new MultidimArray<double>);
	subsetProjections[poolIdx]->resizeNoCopy(xdim, xdim);
	applyGeometry(LINEAR,*(subsetProjections[poolIdx]),Paux(),A,IS_INV,DONT_WRAP,0.);

#ifdef DEBUG
	std::cout << "Row: " << " rot: " << rot << " tilt: " << tilt
			  << " psi: " << psi << " sx: " << x << " sy: " << y
			  << " flip: " << flip << std::endl;
	Image<double> save;
	save()=*(subsetProjections[poolIdx]);
	save.write(formatString("PPPtheo%02d.xmp",poolIdx));
#endif
}

// Select the subset associated to a particleId
void ProgClassifySignificant::selectSubset(size_t particleId, bool &flagEmpty)
{
	size_t poolIdx=0;
	FileName fnImg;
	for (size_t i=0; i<projector.size(); i++)
	{
		subsetAngles[i].clear();
		subsetProjectionIdx[i].clear();
		size_t crIdx=currentRowIdx[i];
		if (crIdx>=setAngles[i].size())
			return;
		MDRow & currentRow=setAngles[i][crIdx];
		/*if (i==0) // First time we see this image
		{
			currentRow.getValue(MDL_IMAGE,fnImg);
			Iexp.read(fnImg);
			std::cout << "Particle fnImg: " << fnImg << std::endl;
		}*/
		size_t currentParticleId;
		currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
		size_t idxMax=setAngles[i].size();
		while (currentParticleId<=particleId)
		{
			if (currentParticleId==particleId)
			{
				flagEmpty=false;
				subsetAngles[i].push_back(currentRow);
				subsetProjectionIdx[i].push_back(poolIdx);
				currentRow.getValue(MDL_IMAGE,fnImg);
				Iexp.push_back(new Image<double>);
				Iexp[poolIdx]->read(fnImg);
				//std::cout << "Particle fnImg: " << fnImg << " in " << poolIdx << std::endl;
				generateProjection(i,poolIdx,currentRow);
				poolIdx+=1;
			}
			crIdx+=1;
			if (crIdx<idxMax)
			{
				currentRow=setAngles[i][crIdx];
				currentRow.getValue(MDL_PARTICLE_ID,currentParticleId);
			}
			else
				break;
		}
		currentRowIdx[i]=crIdx;
	}
#ifdef DEBUG
	std::cout << "Reading " << fnImg << std::endl;
	char c;
	std::cout << "Press any key" << std::endl;
	std::cin >> c;
#endif
}
#undef DEBUG

void computeWeightedCorrelation(MultidimArray<double> &I1, MultidimArray<double> &I2, MultidimArray<double> &Iexp1,
		MultidimArray<double> &Iexp2, double &corr1exp, double &corr2exp, bool I1isEmpty, bool I2isEmpty, int xdim,
		bool onlyIntersection, int numVotes, size_t id, std::ofstream *fs)
{

	MultidimArray<double> Idiff, I2Aligned, Iexp2Aligned;
	Matrix2D<double> M;

	MultidimArray<double> Iaux1(1, 1, xdim, xdim);
	MultidimArray<double> Iaux2(1, 1, xdim, xdim);

	I1.setXmippOrigin();
	I2.setXmippOrigin();
	Iexp1.setXmippOrigin();
	Iexp2.setXmippOrigin();

	I2Aligned=I2;
	Iexp2Aligned=Iexp2;

	if (!I1isEmpty && !I2isEmpty){
		alignImages(I1, I2Aligned, M, false);
		Idiff=I1;
		Idiff-=I2Aligned;
	}
	else if (!I1isEmpty){
		Idiff= I1;
	}
	else if (!I2isEmpty){
		Idiff= -I2;
	}
	/*Image<double> save;
	save()=I1;
	save.write("I1.xmp");
	save()=I2;
	save.write("I2.xmp");
	save()=I2Aligned;
	save.write("I2Aligned.xmp");
	save()=Iexp1;
	save.write("Iexp1.xmp");
	save()=Iexp2;
	save.write("Iexp2.xmp");*/

	double mean, std;
	Idiff.computeAvgStdev(mean,std);
	Idiff.selfABS();
	double threshold=std;

	/*save()=Idiff;
	save.write("Idiff.xmp");*/

	double mean1, std1, mean2, std2, th1, th2;
	I1.computeAvgStdev(mean1,std1);
	th1 = std1;
	I2.computeAvgStdev(mean2,std2);
	th2 = std2;

	//std::cout << "threshold: " << threshold << std::endl;

	corr1exp=corr2exp=0.0;

	if (!I1isEmpty && !I2isEmpty){
		applyGeometry(LINEAR, Iexp2Aligned, Iexp2, M, IS_NOT_INV, false);
	}
	/*save()=Iexp2Aligned;
	save.write("Iexp2Aligned.xmp");*/

	// Estimate the mean and stddev within the mask
	double N=0;
	double N1=0, N2=0;
	double sumMI1=0, sumMI2=0, sumMIexp1=0, sumMIexp2=0;
	double sumI1=0, sumI2=0, sumIexp1=0, sumIexp2=0;
	double sumWI1=0, sumWI2=0, sumWIexp1=0, sumWIexp2=0;

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		double p1=DIRECT_MULTIDIM_ELEM(I1,n);
		double p2=DIRECT_MULTIDIM_ELEM(I2,n);
		double pexp1=DIRECT_MULTIDIM_ELEM(Iexp1,n);
		double pexp2=DIRECT_MULTIDIM_ELEM(Iexp2,n);
		sumI1+=p1;
		sumI2+=p2;
		sumIexp1+=pexp1;
		sumIexp2+=pexp2;
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
		{
			double p2Alg=DIRECT_MULTIDIM_ELEM(I2Aligned,n);
			double pexp2Alg=DIRECT_MULTIDIM_ELEM(Iexp2Aligned,n);
			sumWI1+=p1;
			sumWI2+=p2Alg;
			sumWIexp1+=pexp1;
			sumWIexp2+=pexp2Alg;
			N+=1.0;
			//DIRECT_MULTIDIM_ELEM(Iaux1,n)=1.0;
		}
		if(p1>th1){
			sumMI1+=p1;
			sumMIexp1+=pexp1;
			N1+=1.0;
		}
		if(p2>th2){
			sumMI2+=p2;
			sumMIexp2+=pexp2;
			N2+=1.0;
		}
	}

	//Image<double> save;
	/*save()=Iaux1;
	save.write("Iaux1.xmp");
	save()=Iaux2;
	save.write("Iaux2.xmp");*/


	double sumMI1exp1=0.0, sumMI2exp2=0.0, sumMI1I1=0.0, sumMI2I2=0.0, sumMIexpIexp1=0.0, sumMIexpIexp2=0.0;
	double sumI1exp1=0.0,  sumI2exp2=0.0,  sumI1I1=0.0,  sumI2I2=0.0,  sumIexpIexp1=0.0, sumIexpIexp2=0.0;
	double sumWI1exp1=0.0, sumWI2exp2=0.0, sumWI1I1=0.0, sumWI2I2=0.0, sumWIexpIexp1=0.0, sumWIexpIexp2=0.0;
	double corrM1exp, corrM2exp, corrN1exp, corrN2exp, corrW1exp, corrW2exp;

	double sumWI1exp2=0.0, sumWI2exp1=0.0, corrWI1exp2, corrWI2exp1;
	double sumI1exp2=0.0, sumI2exp1=0.0, corrI1exp2, corrI2exp1;

	double avg1, avgExp1, avgM1, avgMExp1, avgW1, avgWExp1, avg2, avgExp2, avgM2, avgMExp2, avgW2, avgWExp2, iN1, iN2, iN;
	double isize=1.0/MULTIDIM_SIZE(Idiff);
	avg1=sumI1*isize;
	avgExp1=sumIexp1*isize;
	avg2=sumI2*isize;
	avgExp2=sumIexp2*isize;
	if (N1>0){
		iN1=1.0/N1;
		avgM1=sumMI1*iN1;
		avgMExp1=sumMIexp1*iN1;
	}
	if (N2>0){
		iN2=1.0/N2;
		avgM2=sumMI2*iN2;
		avgMExp2=sumMIexp2*iN2;
	}
	if(N>0){
		iN=1.0/N;
		avgW1=sumWI1*iN;
		avgWExp1=sumWIexp1*iN;
		avgW2=sumWI2*iN;
		avgWExp2=sumWIexp2*iN;
	}
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Idiff)
	{
		double p1=DIRECT_MULTIDIM_ELEM(I1,n);
		double p2=DIRECT_MULTIDIM_ELEM(I2,n);
		double pexp1=DIRECT_MULTIDIM_ELEM(Iexp1,n);
		double pexp2=DIRECT_MULTIDIM_ELEM(Iexp2,n);
		double p1a=p1-avg1;
		double p2a=p2-avg2;
		double pexpa1=pexp1-avgExp1;
		double pexpa2=pexp2-avgExp2;
		sumI1exp1+=p1a*pexpa1;
		sumI2exp2+=p2a*pexpa2;
		sumI1I1 +=p1a*p1a;
		sumI2I2 +=p2a*p2a;
		sumIexpIexp1 +=pexp1*pexp1;
		sumIexpIexp2 +=pexp2*pexp2;

		sumI1exp2+=p1a*pexpa2;
		sumI2exp1+=p2a*pexp1;

		if (p1>th1){
			p1a=p1-avgM1;
			pexpa1=pexp1-avgMExp1;
			sumMI1exp1+=p1a*pexpa1;
			sumMI1I1 +=p1a*p1a;
			sumMIexpIexp1 +=pexpa1*pexpa1;
		}
		if (p2>th2){
			p2a=p2-avgM2;
			pexpa2=pexp2-avgMExp2;
			sumMI2exp2+=p2a*pexpa2;
			sumMI2I2 +=p2a*p2a;
			sumMIexpIexp2 +=pexpa2*pexpa2;
		}
		if (DIRECT_MULTIDIM_ELEM(Idiff,n)>threshold)
		{
			double p2Alg=DIRECT_MULTIDIM_ELEM(I2Aligned,n);
			double pexp2Alg=DIRECT_MULTIDIM_ELEM(Iexp2Aligned,n);
			p1a=p1-avgW1;
			p2a=p2Alg-avgW2;
			pexpa1=pexp1-avgWExp1;
			pexpa2=pexp2Alg-avgWExp2;

			double w=DIRECT_MULTIDIM_ELEM(Idiff,n);
			double wp1a=w*p1a;
			double wp2a=w*p2a;
			double wpexpa1=w*pexpa1;
			double wpexpa2=w*pexpa2;

			sumWI1exp1+=wp1a*pexpa1;
			sumWI2exp2+=wp2a*pexpa2;
			sumWI1I1 +=wp1a*p1a;
			sumWI2I2 +=wp2a*p2a;
			sumWIexpIexp1 +=wpexpa1*pexpa1;
			sumWIexpIexp2 +=wpexpa2*pexpa2;

			sumWI1exp2+=wp1a*pexpa2;
			sumWI2exp1+=wp2a*pexpa1;
		}
	}

	sumI1exp1*=isize;
	sumI2exp2*=isize;
	sumI1I1*=isize;
	sumI2I2*=isize;
	sumIexpIexp1*=isize;
	sumIexpIexp2*=isize;

	sumI1exp2*=isize;
	sumI2exp1*=isize;

	sumMI1exp1*=iN1;
	sumMI2exp2*=iN2;
	sumMI1I1*=iN1;
	sumMI2I2*=iN2;
	sumMIexpIexp1*=iN1;
	sumMIexpIexp2*=iN2;

	sumWI1exp1*=iN;
	sumWI2exp2*=iN;
	sumWI1I1*=iN;
	sumWI2I2*=iN;
	sumWIexpIexp1*=iN;
	sumWIexpIexp2*=iN;

	sumWI1exp2*=iN;
	sumWI2exp1*=iN;

	corrM1exp=sumMI1exp1/sqrt(sumMI1I1*sumMIexpIexp1);
	corrM2exp=sumMI2exp2/sqrt(sumMI2I2*sumMIexpIexp2);
	corrN1exp=sumI1exp1/sqrt(sumI1I1*sumIexpIexp1);
	corrN2exp=sumI2exp2/sqrt(sumI2I2*sumIexpIexp2);
	corrW1exp=sumWI1exp1/sqrt(sumWI1I1*sumWIexpIexp1);
	corrW2exp=sumWI2exp2/sqrt(sumWI2I2*sumWIexpIexp2);

	corrWI1exp2 = sumWI1exp2/sqrt(sumWI1I1*sumWIexpIexp2);
	corrWI2exp1 = sumWI2exp1/sqrt(sumWI2I2*sumWIexpIexp1);

	corrI1exp2 = sumI1exp2/sqrt(sumI1I1*sumIexpIexp2);
	corrI2exp1 = sumI2exp1/sqrt(sumI2I2*sumIexpIexp1);

	if (isnan(corrWI2exp1))
		corrWI2exp1=-1.0;
	if (isnan(corrWI1exp2))
		corrWI1exp2=-1.0;

	if(isnan(corrN1exp))
		corrN1exp=-1.0;
	if(isnan(corrM1exp))
		corrM1exp=-1.0;
	if(isnan(corrW1exp))
		corrW1exp=-1.0;
	if(isnan(corrN2exp))
		corrN2exp=-1.0;
	if(isnan(corrM2exp))
		corrM2exp=-1.0;
	if(isnan(corrW2exp))
		corrW2exp=-1.0;
	if (isnan(corrI2exp1))
		corrI2exp1=-1.0;
	if (isnan(corrI1exp2))
		corrI1exp2=-1.0;
	corr2exp = -1;
	corr1exp=-1;
	if (onlyIntersection){
		if (corrN1exp==-1 || corrN2exp==-1)
			return;
	}
	double imedN1exp=imedDistance(I1, Iexp1);
	double imedN2exp=imedDistance(I2, Iexp2);

	int votes=0;
	if (corrN1exp>corrN2exp && corrN1exp>0)
		votes+=1;
	else if (corrN1exp<corrN2exp && corrN2exp>0)
		votes-=1;
	if (corrM1exp>corrM2exp && corrM1exp>0)
		votes+=1;
	else if (corrM1exp<corrM2exp && corrM2exp>0)
		votes-=1;
	if (corrW1exp>corrW2exp && corrW1exp>0)
		votes+=1;
	else if (corrW1exp<corrW2exp && corrW2exp>0)
		votes-=1;
	if ((corrW1exp-corrWI2exp1)>(corrW2exp-corrWI1exp2))
		votes+=1;
	else if ((corrW1exp-corrWI2exp1)<(corrW2exp-corrWI1exp2))
		votes-=1;
	if (imedN1exp<imedN2exp)
		votes+=1;
	else if (imedN1exp>imedN2exp)
		votes-=1;
	if (corrN1exp<corrI2exp1)
		votes-=1;
	else if (corrN2exp<corrI1exp2)
		votes+=1;
//	if ((corrN1exp-corrI2exp1)>(corrN2exp-corrI1exp2))
//		votes+=1;
//	else if ((corrN1exp-corrI2exp1)<(corrN2exp-corrI1exp2))
//		votes-=1;

	if(votes>=numVotes){
		corr1exp = corrN1exp;
		corr2exp = -1;
	}
	else if(votes<=-numVotes){
		corr2exp= corrN2exp;
		corr1exp=-1;
	}

	//(*fs) << id << " " << corr1exp << " " << corrN1exp << " " << corrM1exp << " " << corrW1exp << " " << corrI2exp1 << " " << corrWI2exp1 << " " << corr2exp << " " << corrN2exp << " " << corrM2exp << " " << corrW2exp << " " << corrI1exp2 << " " << corrWI1exp2 << std::endl;
//	std::cout << "corr1exp= " << corr1exp << " corrN1exp: " << corrN1exp << " corrM1exp=" << corrM1exp << " corrW1exp=" << corrW1exp << " corrWI2exp1=" << corrWI2exp1 << " corrI2exp1=" << corrI2exp1 << " imedN1exp=" << imedN1exp << std::endl;
//	std::cout << "corr2exp= " << corr2exp << " corrN2exp: " << corrN2exp << " corrM2exp=" << corrM2exp << " corrW2exp=" << corrW2exp << " corrWI1exp2=" << corrWI1exp2 << " corrI1exp2=" << corrI1exp2 << " imedN2exp=" << imedN2exp << std::endl;
//	std::cout << "votes= " << votes << std::endl;
}

void ProgClassifySignificant::updateClass(int n, double wn)
{
	double CCbest=-1e38;
	int iCCbest=-1;
	VMetaData &subsetAngles_n=subsetAngles[n];
	for (int i=0; i<subsetAngles_n.size(); i++)
	{
		double cc;
		subsetAngles_n[i].getValue(MDL_MAXCC,cc);
		if (cc>CCbest)
		{
			CCbest=cc;
			iCCbest=i;
		}
	}
	if (iCCbest>=0)
	{
		MDRow newRow=subsetAngles_n[iCCbest];
		// COSS newRow.setValue(MDL_WEIGHT,wn);
		classifiedAngles[n].push_back(newRow);
	}
}

//#define DEBUG
void ProgClassifySignificant::run()
{
	show();
	produceSideInfo();

	std::ofstream fs("./correlations.txt", std::ofstream::out);

	if (verbose>0)
	{
		std::cerr << "Classifying images ..." << std::endl;
		init_progress_bar(setIds.size());
	}

	Matrix1D<double> winning(projector.size());
	Matrix1D<double> corrDiff(projector.size());
	Matrix1D<double> weight;

	MultidimArray<double> I1, I2;
	MultidimArray<double> Iexp1, Iexp2;
	bool I1isEmpty, I2isEmpty;
	double corr1exp, corr2exp;
	bool flagEmpty=true;

	for (size_t iid=0; iid<setIds.size(); iid++)
	{
		flagEmpty=true;
		size_t id=setIds[iid];
		selectSubset(id, flagEmpty);
		if (flagEmpty)
			continue;

		winning.initZeros();
		corrDiff.initZeros();
		Image<double> Iaux = *Iexp[0];
		int xdim = (int)XSIZE(Iaux());
		I1.initZeros(1, 1, xdim, xdim);
		I2.initZeros(1, 1, xdim, xdim);
		Iexp1.initZeros(1, 1, xdim, xdim);
		Iexp2.initZeros(1, 1, xdim, xdim);

		for (size_t ivol1=0; ivol1<projector.size(); ivol1++)
		{
			std::vector<size_t> &subset1=subsetProjectionIdx[ivol1];
			for (size_t ivol2=ivol1+1; ivol2<projector.size(); ivol2++)
			{
				std::vector<size_t> &subset2=subsetProjectionIdx[ivol2];
				size_t i1=0;
				do //for (size_t i1=0; i1<subset1.size(); i1++)
				{
					//std::cout << "subset1.size()" << subset1.size() << std::endl;
					if (subset1.size()==0)
						I1isEmpty=true;
					else
					{
						//MultidimArray<double> &I1=*(subsetProjections[subset1[i1]]);
						I1=*(subsetProjections[subset1[i1]]);
						Iexp1 = (*(Iexp[subset1[i1]]))();
						I1isEmpty = false;
						//std::cout << "subset1[i1]" << subset1[i1] << std::endl;
					}

					size_t i2=0;
					do //for (size_t i2=0; i2<subset2.size(); i2++)
					{
						//std::cout << "subset2.size()" << subset2.size() << std::endl;
						if (subset2.size()==0)
							I2isEmpty=true;
						else
						{
							//MultidimArray<double> &I2=*(subsetProjections[subset2[i2]]);
							I2=*(subsetProjections[subset2[i2]]);
							Iexp2 = (*(Iexp[subset2[i2]]))();
							I2isEmpty=false;
							//std::cout << "subset2[i2]" << subset2[i2] << std::endl;
						}

						//if (!I1isEmpty && !I2isEmpty)
						//{
							computeWeightedCorrelation(I1, I2, Iexp1, Iexp2, corr1exp, corr2exp, I1isEmpty, I2isEmpty,
									xdim, onlyIntersection, numVotes, id, &fs);
							if (isnan(corr1exp))
								corr1exp=-1.0;
							if (isnan(corr2exp))
								corr2exp=-1.0;

							double corrDiff12=corr1exp-corr2exp;
							if (corrDiff12>0 && corr1exp>0)
							{
								VEC_ELEM(winning,ivol1)+=1;
								VEC_ELEM(corrDiff,ivol1)+=corrDiff12;
								VEC_ELEM(corrDiff,ivol2)-=corrDiff12;
							}
							else if (corrDiff12<0 && corr2exp>0)
							{
								VEC_ELEM(winning,ivol2)+=1;
								VEC_ELEM(corrDiff,ivol2)-=corrDiff12;
								VEC_ELEM(corrDiff,ivol1)+=corrDiff12;
							}
						//}
						/*else if (!I1isEmpty)
						{

							VEC_ELEM(winning,ivol1)+=1;
							VEC_ELEM(corrDiff,ivol1)+=corr1exp;
							VEC_ELEM(corrDiff,ivol2)-=corr1exp;
						}
						else if (!I2isEmpty)
						{
							VEC_ELEM(winning,ivol2)+=1;
							VEC_ELEM(corrDiff,ivol1)-=corr2exp;
							VEC_ELEM(corrDiff,ivol2)+=corr2exp;
						}*/

						/*Image<double> save;
						save()=Iexp1;
						save.write("PPPexp1.xmp");
						save()=Iexp2;
						save.write("PPPexp2.xmp");
						save()=I1;
						save.write("PPP1.xmp");
						save()=I2;
						save.write("PPP2.xmp");
						std::cout << id << " " << corr1exp << " " << corr2exp << " " << I1isEmpty << " " << I2isEmpty << std::endl;
						std::cout << "winning=" << winning << std::endl;
						std::cout << "corrDiff=" << corrDiff << std::endl;

						char c;
						std::cout << "Press any key" << std::endl;
						std::cin >> c;*/

						i2++;
					}while(i2<subset2.size());

					i1++;
				}while(i1<subset1.size());
			}
		}
		double iNcomparisons=1.0/winning.sum();
		winning*=iNcomparisons;
		weight=corrDiff;
		weight*=iNcomparisons;
		weight*=winning;

//		std::cout << corrDiff << std::endl;
//		std::cout << winning << std::endl;
//		std::cout << weight << std::endl;


		int nBest=weight.maxIndex();
		double wBest=VEC_ELEM(weight,nBest);
		if (wBest>0)
			updateClass(nBest,wBest);
		if (verbose>0)
			progress_bar(iid);

	}
	progress_bar(setIds.size());

	// Write output
	MetaData md;
	for (size_t ivol=0; ivol<projector.size(); ivol++)
	{
		size_t objId=md.addObject();
		md.setValue(MDL_REF3D,(int)ivol+1,objId);
		md.setValue(MDL_CLASS_COUNT,classifiedAngles[ivol].size(),objId);
	}
	md.write("classes@"+fnOut);
	for (size_t ivol=0; ivol<projector.size(); ivol++)
	{
		md.clear();
		if (classifiedAngles[ivol].size()>0)
		{
			md.fromVMetaData(classifiedAngles[ivol]);
			double currentWmax=md.getColumnMax(MDL_WEIGHT);
			double currentWmin=md.getColumnMin(MDL_WEIGHT);
			if (currentWmax>currentWmin)
				md.operate(formatString("weight=%f*(weight-%f)+%f",(1.0-wmin)/(currentWmax-currentWmin),currentWmin,wmin));
			else
				md.operate(formatString("weight=%f",wmin));
			md.setValueCol(MDL_REF3D,(int)ivol+1);
		}
		else
			REPORT_ERROR(ERR_VALUE_EMPTY,formatString("Class %d have been depleted of images. Cannot continue processing",ivol));
		md.write(formatString("class%06d_images@%s",ivol+1,fnOut.c_str()),MD_APPEND);
	}

	fs.close();

}

#undef DEBUG
