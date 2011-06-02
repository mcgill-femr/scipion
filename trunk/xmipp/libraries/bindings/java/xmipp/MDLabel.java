package xmipp; 
public class MDLabel {
   public static final int MDL_UNDEFINED = -1;
   public static final int MDL_FIRST_LABEL = 0;  ///< The label MDL_OBJID is special and should not be used
   public static final int MDL_OBJID = MDL_FIRST_LABEL; ///< object id (int); NOTE: This label is special and shouldn't be used
   public static final int MDL_ANGLE_COMPARISON = 1;  ///< Angular comparison (see angular_distance.cpp)
   public static final int MDL_ANGLEPSI2 = 2;  ///< Psi angle of an image (double = 2; degrees)
   public static final int MDL_ANGLEPSI = 3;  ///< Psi angle of an image (double = 3; degrees)
   public static final int MDL_ANGLEROT2 = 4;  ///< Rotation angle of an image (double = 4; degrees)
   public static final int MDL_ANGLEROT = 5;  ///< Rotation angle of an image (double = 5; degrees)
   public static final int MDL_ANGLETILT2 = 6;  ///< Tilting angle of an image (double = 6; degrees)
   public static final int MDL_ANGLETILT = 7;  ///< Tilting angle of an image (double = 7; degrees)
   public static final int MDL_ASSOCIATED_IMAGE1 = 8;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE2 = 9;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE3 = 10;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE4 = 11;  ///< Image associated to this object (std::string)
   public static final int MDL_ASSOCIATED_IMAGE5 = 12;  ///< Image associated to this object (std::string)
   public static final int MDL_AVG = 13;  ///< average value (double)
   public static final int MDL_AZIMUTALANGLE = 14;  ///< ctf definition azimutal angle
   public static final int MDL_BGMEAN = 15;  ///< Mean background value for an image
   public static final int MDL_BLOCK = 16;  ///< Current block number (for incremental EM)
   public static final int MDL_CELLX = 17;  ///< Cell location for crystals
   public static final int MDL_CELLY = 18;  ///< Cell location for crystals
   public static final int MDL_CLASSIFICATION_DATA = 19;  ///< Data vector for classification (vector double)
   public static final int MDL_CLASSIFICATION_DATA_SIZE = 20;  ///< Size of data vectors for classification (int)
   public static final int MDL_CLASSIFICATION_INTRACLASS_DISTANCE = 21;  ///< Average intraclass distance (double)
   public static final int MDL_COMMENT = 22;  ///< A comment for this object /*** NOTE THIS IS A SPECIAL CASE AND SO IS TREATED ***/
   public static final int MDL_COST = 23;  ///< Cost for the image (double)
   public static final int MDL_COUNT = 24;  ///< Number of elements of a type (int) [this is a genereic type do not use to transfer information to another program]
   public static final int MDL_CTFINPUTPARAMS = 25;  ///< Parameters file for the CTF Model (std::string)
   public static final int MDL_CTFMODEL = 26;  ///< Name for the CTF Model (std::string)
   public static final int MDL_CTFMODEL2 = 27;  ///< Name for another CTF model (std::string)
   public static final int MDL_CTF_SAMPLING_RATE = 28;  ///< Sampling rate
   public static final int MDL_CTF_SAMPLING_RATE_Z = 29;  ///< Sampling rate in Z direction
   public static final int MDL_CTF_VOLTAGE = 30;  ///< Microscope voltage (kV)
   public static final int MDL_CTF_DEFOCUSA = 31;  ///< aver (Angage defocusstroms)
   public static final int MDL_CTF_DEFOCUSU = 32;  ///< Defocus U (Angstroms)
   public static final int MDL_CTF_DEFOCUSV = 33;  ///< Defocus V (Angstroms)
   public static final int MDL_CTF_DEFOCUS_ANGLE = 34;  ///< Defocus angle (degrees)
   public static final int MDL_CTF_CS = 35;  ///< Spherical aberration
   public static final int MDL_CTF_CA = 36;  ///< Chromatic aberration
   public static final int MDL_CTF_GROUP = 37;  ///< group images by defocus
   public static final int MDL_CTF_ENERGY_LOSS = 38;  ///< Energy loss
   public static final int MDL_CTF_LENS_STABILITY = 39;  ///< Lens stability
   public static final int MDL_CTF_CONVERGENCE_CONE = 40;  ///< Convergence cone
   public static final int MDL_CTF_LONGITUDINAL_DISPLACEMENT = 41;  ///< Longitudinal displacement
   public static final int MDL_CTF_TRANSVERSAL_DISPLACEMENT = 42;  ///< Transversal displacemente
   public static final int MDL_CTF_Q0 = 43;  ///< Inelastic absorption
   public static final int MDL_CTF_K = 44;  ///< CTF gain
   public static final int MDL_CTFBG_GAUSSIAN_K = 45;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_SIGMAU = 46;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_SIGMAV = 47;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_CU = 48;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_CV = 49;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN_ANGLE = 50;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_K = 51;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_U = 52;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_V = 53;  ///< CTF Background parameter
   public static final int MDL_CTFBG_SQRT_ANGLE = 54;  ///< CTF Background parameter
   public static final int MDL_CTFBG_BASELINE = 55;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_K = 56;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_SIGMAU = 57;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_SIGMAV = 58;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_CU = 59;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_CV = 60;  ///< CTF Background parameter
   public static final int MDL_CTFBG_GAUSSIAN2_ANGLE = 61;  ///< CTF Background parameter
   public static final int MDL_CTF_CRITERION_PSDCORRELATION90 = 62;  ///< PSD correlation at 90 degrees
   public static final int MDL_CTF_CRITERION_FIRSTZERORATIO = 63;  ///< First zero ratio
   public static final int MDL_CTF_CRITERION_FIRSTZEROAVG = 64;  ///< First zero average (in Angstroms)
   public static final int MDL_CTF_CRITERION_FIRSTZERODISAGREEMENT = 65;  ///< First zero disagreement with second model (in Angstroms)
   public static final int MDL_CTF_CRITERION_DAMPING = 66;  ///< Minimum damping at border
   public static final int MDL_CTF_CRITERION_PSDRADIALINTEGRAL = 67;  ///< Integral of the radial PSD
   public static final int MDL_CTF_CRITERION_FITTINGSCORE = 68;  ///< Score of the fitting
   public static final int MDL_CTF_CRITERION_FITTINGCORR13 = 69;  ///< Correlation between the 1st and 3rd ring of the CTF
   public static final int MDL_CTF_CRITERION_PSDVARIANCE = 70;  ///< PSD variance
   public static final int MDL_CTF_CRITERION_PSDPCA1VARIANCE = 71;  ///< Variance in the first principal component of the PSDs
   public static final int MDL_CTF_CRITERION_PSDPCARUNSTEST = 72;  ///< Runs test on the projection of the PSD on the first principal component
   public static final int MDL_CTF_CRITERION_COMBINED = 73;  ///< Combined criterion formed by several other criteria
   public static final int MDL_CTF_CRITERION_NORMALITY = 74;  ///< Normality test between histogram of micrography and gaussian distribution
   public static final int MDL_CTF_XRAY_DIMENSIONS = 75;  // Size in pixels of the 3D PSF to be created (Xdim = 75;  Ydim = 75;  Zdim)
   public static final int MDL_CTF_XRAY_LAMBDA = 76;  /// X-ray wavelength (nm)
   public static final int MDL_CTF_XRAY_LENS_TYPE = 77;  ///Algorithm used to generate Xray PSF
   public static final int MDL_CTF_XRAY_MAGNIFICATION = 78;  /// Magnification of the X-ray microscope
   public static final int MDL_CTF_XRAY_OUTER_ZONE_WIDTH = 79;  /// Outermost zone width of the X-ray Fresnel lens (nm)
   public static final int MDL_CTF_XRAY_ZONES_NUMBER = 80;  // Number of zones of the X-ray Fresnel lens
   public static final int MDL_DATATYPE = 81;  ///< if read from file original image datatype = 81;  this is an struct defined in image
   public static final int MDL_DEFGROUP = 82;  ///< Defocus group
   public static final int MDL_DIRECTION = 83;  ///< Direction in 3D
   public static final int MDL_DM3_IDTAG = 84; 
   public static final int MDL_DM3_NODEID = 85; 
   public static final int MDL_DM3_NUMBER_TYPE = 86; 
   public static final int MDL_DM3_PARENTID = 87; 
   public static final int MDL_DM3_TAGCLASS = 88; 
   public static final int MDL_DM3_TAGNAME = 89; 
   public static final int MDL_DM3_SIZE = 90; 
   public static final int MDL_DM3_VALUE = 91; 
   public static final int MDL_ENABLED = 92;  ///< Is this image enabled? (int [-1 or 1])
   public static final int MDL_FLIP = 93;  ///< Flip the image? (bool)
   public static final int MDL_IMAGE_CLASS_COUNT = 94;  ///< Number of images assigned to the same class as this image
   public static final int MDL_IMAGE_CLASS_GROUP = 95;  ///< Name of the class group for this image (metadata with all the images assigned to that class)
   public static final int MDL_IMAGE_CLASS = 96;  ///< Name of the class representative for this image
   public static final int MDL_IMAGE = 97;  ///< Name of an image (std::string)
   public static final int MDL_IMAGE_ORIGINAL = 98;  ///< Name of an image from which MDL_IMAGE is coming from
   public static final int MDL_IMAGE_TILTED = 99;  ///< Name of the tilted images associated to MDL_IMAGE
   public static final int MDL_IMGMD = 100;  ///< Name of Metadata file for all images (string)
   public static final int MDL_INTSCALE = 101;  ///< Intensity scale for an image
   public static final int MDL_ITER = 102;  ///< Current iteration number (int)
   public static final int MDL_K = 103;  ///< //ctf definition K
   public static final int MDL_KERDENSOM_FUNCTIONAL = 104;  ///< Functional value (double)
   public static final int MDL_KERDENSOM_REGULARIZATION = 105;  ///< Regularization value (double)
   public static final int MDL_KERDENSOM_SIGMA = 106;  ///< Sigma value (double)
   public static final int MDL_KSTEST = 107;  ///<KS-test statistics
   public static final int MDL_LL = 108;  ///< contribution of an image to log-likelihood value
   public static final int MDL_MAPTOPOLOGY = 109;  ///< Map topology (KerDenSOM = 109;  ...)
   public static final int MDL_MASK = 110;  ///< Name of a mask associated to image
   public static final int MDL_MAXCC = 111;  ///< Cross-correlation for the image (double)
   public static final int MDL_MAX = 112;  ///<maximum value (double)
   public static final int MDL_MICROGRAPH = 113;  ///< Name of a micrograph (std::string)
   public static final int MDL_MICROGRAPH_TILTED = 114;  ///< Name of the corresponding tilted micrograph (std::string)
   public static final int MDL_MIN = 115;  ///<minimum value (double)
   public static final int MDL_MIRRORFRAC = 116;  ///< Mirror fraction for a Maximum Likelihood model
   public static final int MDL_MISSINGREGION_NR = 117;  ///< Number of missing region in subtomogram
   public static final int MDL_MISSINGREGION_TYPE = 118;  ///< Type of missing region in subtomogram
   public static final int MDL_MISSINGREGION_THY0 = 119;  ///< Initial tilt angle in Y for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THYF = 120;  ///< Final tilt angle in Y for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THX0 = 121;  ///< Initial tilt angle in X for missing region in subtomogram
   public static final int MDL_MISSINGREGION_THXF = 122;  ///< Final tilt angle in X for missing region in subtomogram
   public static final int MDL_MODELFRAC = 123;  ///< Model fraction (alpha_k) for a Maximum Likelihood model
   public static final int MDL_NEIGHBORS = 124;  ///< Vector of indexes to points some "neighbors"
   public static final int MDL_NMA = 125;  ///< Normal mode displacements (vector double)
   public static final int MDL_NMA_MODEFILE = 126;  ///< File with an NMA mode
   public static final int MDL_NOISE_ANGLES = 127;  ///< Noise description for projected angles
   public static final int MDL_NOISE_PARTICLE_COORD = 128;  ///< Noise description for particle's center coordenates (when projecting)
   public static final int MDL_NOISE_PIXEL_LEVEL = 129;  ///< Noise description for pixels' gray level (when projecting)
   public static final int MDL_ORDER = 130;  /// auxiliary label to be used as an index (long)
   public static final int MDL_ORIGINX = 131;  ///< Origin for the image in the X axis (double)
   public static final int MDL_ORIGINY = 132;  ///< Origin for the image in the Y axis (double)
   public static final int MDL_ORIGINZ = 133;  ///< Origin for the image in the Z axis (double)
   public static final int MDL_PMAX = 134;  ///< Maximum value of normalized probability function (now called "Pmax/sumP") (double)
   public static final int MDL_PRJ_DIMENSIONS = 135;  // X = 135; Y dimensions for the generated projections
   public static final int MDL_PRJ_TILT_RANGE = 136;  // Vector with the initial and final tilt angle values = 136;  and step size
   public static final int MDL_PRJ_VOL = 137;         // Volume file name to generate projections from
   public static final int MDL_PSD = 138;  ///< A Power Spectrum Density file name (std::string)
   public static final int MDL_RANDOMSEED = 139;  ///< Seed for random number generator
   public static final int MDL_REF3D = 140;  ///< 3D Class to which the image belongs (int)
   public static final int MDL_REF = 141;  ///< Class to which the image belongs (int)
   public static final int MDL_REFMD = 142;  ///< Name of Metadata file for all references(string)
   public static final int MDL_RESOLUTION_DPR = 143;  ///<differential phase residual (double)
   public static final int MDL_RESOLUTION_ERRORL2 = 144;  ///<Error in l2 (double)
   public static final int MDL_RESOLUTION_FRC = 145;  ///<Fourier shell correlation (double)
   public static final int MDL_RESOLUTION_FRCRANDOMNOISE = 146;  ///<Fourier shell correlation noise (double)
   public static final int MDL_RESOLUTION_FREQ = 147;  ///<Frequency in 1/A (double)
   public static final int MDL_RESOLUTION_FREQREAL = 148;  ///< Frequency in A (double)
   public static final int MDL_SAMPLINGRATE = 149;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEX = 150;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEY = 151;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SAMPLINGRATEZ = 152;  ///< sampling rate in A/pixel (double)
   public static final int MDL_SCALE = 153;  ///< scaling factor for an image or volume (double)
   public static final int MDL_SELFILE = 154;  ///< Name of an image (std::string)
   public static final int MDL_SERIE = 155;  ///< A collection of micrographs = 155;  e.g. a tilt serie (std::string)
   public static final int MDL_SHIFTX = 156;  ///< Shift for the image in the X axis (double)
   public static final int MDL_SHIFTY = 157;  ///< Shift for the image in the Y axis (double)
   public static final int MDL_SHIFTZ = 158;  ///< Shift for the image in the Z axis (double)
   public static final int MDL_SHIFT_CRYSTALX = 159;  ///< Shift for the image in the X axis (double) for crystals
   public static final int MDL_SHIFT_CRYSTALY = 160;  ///< Shift for the image in the Y axis (double) for crystals
   public static final int MDL_SHIFT_CRYSTALZ = 161;  ///< Shift for the image in the Z axis (double) for crystals
   public static final int MDL_SIGMANOISE = 162;  ///< Standard deviation of the noise in ML model
   public static final int MDL_SIGMAOFFSET = 163;  ///< Standard deviation of the offsets in ML model
   public static final int MDL_SIGNALCHANGE = 164;  ///< Signal change for an image
   public static final int MDL_SPHERICALABERRATION = 165;  ///<ctf definition azimutal angle
   public static final int MDL_STDDEV = 166;  ///<stdandard deviation value (double)
   public static final int MDL_SUM = 167;  ///< Sum of elements of a given type (double) [this is a genereic type do not use to transfer information to another program]
   public static final int MDL_SUMWEIGHT = 168;  ///< Sum of all weights in ML model
   public static final int MDL_SYMNO = 169;  ///< Symmetry number for a projection (used in ART)
   public static final int MDL_TRANSFORMATIONMTRIX = 170;  ///< transformation matrix(vector double)
   public static final int MDL_VOLTAGE = 171;  ///< microscope voltage (double)
   public static final int MDL_WEIGHT = 172;  ///< Weight assigned to the image (double)
   public static final int MDL_WROBUST = 173;  ///< Weight of t-student distribution in robust Maximum likelihood
   public static final int MDL_X = 174;  ///< X component (double)
   public static final int MDL_XINT = 175;  ///< X component (int)
   public static final int MDL_XINTTILT = 176;  ///< X component in tilted micrograph (int)
   public static final int MDL_XSIZE = 177;  ///< X size (int)
   public static final int MDL_Y = 178;  ///< Y component (double)
   public static final int MDL_YINT = 179;  ///< Y component (int)
   public static final int MDL_YINTTILT = 180;  ///< Y component in tilted micrograph (int)
   public static final int MDL_YSIZE = 181;  ///< Y size (int)
   public static final int MDL_Z = 182;  ///< Z component (double)
   public static final int MDL_ZINT = 183;  ///< Z component (int)
   public static final int MDL_ZSCORE = 184;  ///< Z Score (double)
   public static final int MDL_ZSIZE = 185;  ///< Z size (int)
   public static final int MDL_LAST_LABEL = 186;                       // **** NOTE ****: Do keep this label always at the end
}
