/***************************************************************************
 *
 * Authors:     Javier Angel Velazquez Muriel (javi@cnb.uam.es)
 *              Carlos Oscar S�nchez Sorzano (coss@cnb.uam.es)
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
 *  e-mail address 'xmipp@cnb.uam.es'                                  
 ***************************************************************************/

#include "../Prog_adjust_CTF.hh"
#include <XmippData/xmippArgs.hh>
#include <XmippData/xmippHistograms.hh>

/* prototypes */
double CTF_fitness(double *);

/* Number of CTF parameters */
#define CTF_PARAMETERS             24
#define PARAMETRIC_CTF_PARAMETERS  13
#define BACKGROUND_CTF_PARAMETERS  11
#define SQRT_CTF_PARAMETERS         5
#define ENVELOPE_PARAMETERS         8
#define DEFOCUS_PARAMETERS          3
#define FIRST_SQRT_PARAMETER       13
#define FIRST_ENVELOPE_PARAMETER    4
#define FIRST_DEFOCUS_PARAMETER     0

/* Global variables -------------------------------------------------------- */
// Some aliases
Adjust_CTF_Parameters *global_prm;
matrix2D<double>      *f;               // The CTF to model
matrix1D<double>      *global_adjust;   // Current theoretical adjustment

// Frequency of each point in digital units
matrix2D<double>       global_x_digfreq;
matrix2D<double>       global_y_digfreq;
matrix2D<double>       global_w_digfreq;
matrix2D<double>       global_x_contfreq;
matrix2D<double>       global_y_contfreq;
matrix2D<double>       global_w_contfreq;
matrix2D<double>       global_mask;

// Penalization for forbidden values of the parameters
double		       global_heavy_penalization; 

// Penalization factor for the background
bool                   global_penalize;
double                 global_current_penalty;
const double           global_penalty=32; // Maximum penalization

// Speed up factor
int                    global_evaluation_reduction;

// CTF model and noise model
XmippCTF               global_ctfmodel;

// Maximum of the gaussian
double                 global_max_gauss_freq; 

// Autofocus
double                 global_value_th;
double                 global_min_freq;
double                 global_max_freq;

// Program status
int                    global_action; // 0: Computing the background (sqrt)
                                      // 1: Computing the full background
                                      // 2: Computing the envelope
                                      // 3: Computing defoci
                                      // 4: Computing all parameters
                                      // 5: Produce output
int                    global_show;   // 0: Do not show
                                      // 1: Partially detailed
                                      // 2: Very detailed

/* Assign ctfmodel from a vector and viceversa ----------------------------- */
void assign_CTF_from_parameters(double *p, XmippCTF &ctfmodel,
   int ia, int l, bool astigmatic_noise) {
   ctfmodel.Tm=global_prm->Tm; 
   if (ia<= 0 && l>0) {ctfmodel.DeltafU        =p[ 0]; l--;}
   if (ia<= 1 && l>0) {ctfmodel.DeltafV        =p[ 1]; l--;}
   if (ia<= 2 && l>0) {ctfmodel.azimuthal_angle=p[ 2]; l--;}
   if (ia<= 3 && l>0) {ctfmodel.kV             =p[ 3]; l--;}
   if (ia<= 4 && l>0) {ctfmodel.K              =p[ 4]; l--;}
   if (ia<= 5 && l>0) {ctfmodel.Cs             =p[ 5]; l--;}
   if (ia<= 6 && l>0) {ctfmodel.Ca             =p[ 6]; l--;}
   if (ia<= 7 && l>0) {ctfmodel.espr           =p[ 7]; l--;}
   if (ia<= 8 && l>0) {ctfmodel.ispr           =p[ 8]; l--;}
   if (ia<= 9 && l>0) {ctfmodel.alpha          =p[ 9]; l--;}
   if (ia<=10 && l>0) {ctfmodel.DeltaF         =p[10]; l--;}
   if (ia<=11 && l>0) {ctfmodel.DeltaR         =p[11]; l--;}
   if (ia<=12 && l>0) {ctfmodel.Q0             =p[12]; l--;}
   if (ia<=13 && l>0) {ctfmodel.base_line      =p[13]; l--;}        //     0
   if (ia<=14 && l>0) {ctfmodel.sqrt_K         =p[14]; l--;}        //     1
   if (ia<=15 && l>0) {ctfmodel.sqU            =p[15]; l--;}        //     2
   if (ia<=16 && l>0)
      if (astigmatic_noise) {ctfmodel.sqV           =p[16]; l--;}   //     3 *
      else                  {ctfmodel.sqV           =p[15]; l--;}
   if (ia<=17 && l>0)
      if (astigmatic_noise) {ctfmodel.sqrt_angle    =p[17]; l--;}   //     4 *
      else                  {ctfmodel.sqrt_angle    =0;     l--;}
   if (ia<=18 && l>0) {ctfmodel.gaussian_K     =p[18]; l--;}        //     5
   if (ia<=19 && l>0) {ctfmodel.sigmaU         =p[19]; l--;}        //     6
   if (ia<=20 && l>0)
      if (astigmatic_noise) {ctfmodel.sigmaV=p[20]; l--;}           //     7 *
      else                  {ctfmodel.sigmaV=p[19]; l--;}
   if (ia<=21 && l>0)                                               
      if (astigmatic_noise) {ctfmodel.gaussian_angle=p[21]; l--;}   //     8 *
      else                  {ctfmodel.gaussian_angle=0;     l--;}
   if (ia<=22 && l>0) {ctfmodel.cU             =p[22]; l--;}        //     9
   if (ia<=23 && l>0)
      if (astigmatic_noise) {ctfmodel.cV            =p[23]; l--;}   //    10 *
      else                  {ctfmodel.cV            =p[22]; l--;}
}

void assign_parameters_from_CTF(XmippCTF &ctfmodel, double *p,
   int ia, int l, bool astigmatic_noise) {
   if (ia<= 0 && l>0) {p[ 0]=ctfmodel.DeltafU; l--;}
   if (ia<= 1 && l>0) {p[ 1]=ctfmodel.DeltafV; l--;}
   if (ia<= 2 && l>0) {p[ 2]=ctfmodel.azimuthal_angle; l--;}
   if (ia<= 3 && l>0) {p[ 3]=ctfmodel.kV; l--;}
   if (ia<= 4 && l>0) {p[ 4]=ctfmodel.K; l--;}
   if (ia<= 5 && l>0) {p[ 5]=ctfmodel.Cs; l--;}
   if (ia<= 6 && l>0) {p[ 6]=ctfmodel.Ca; l--;}
   if (ia<= 7 && l>0) {p[ 7]=ctfmodel.espr; l--;}
   if (ia<= 8 && l>0) {p[ 8]=ctfmodel.ispr; l--;}
   if (ia<= 9 && l>0) {p[ 9]=ctfmodel.alpha; l--;}
   if (ia<=10 && l>0) {p[10]=ctfmodel.DeltaF; l--;}
   if (ia<=11 && l>0) {p[11]=ctfmodel.DeltaR; l--;}
   if (ia<=12 && l>0) {p[12]=ctfmodel.Q0; l--;}
   if (ia<=13 && l>0) {p[13]=ctfmodel.base_line; l--;}
   if (ia<=14 && l>0) {p[14]=ctfmodel.sqrt_K; l--;}
   if (ia<=15 && l>0) {p[15]=ctfmodel.sqU; l--;}
   if (ia<=16 && l>0)
      if (astigmatic_noise) {p[16]=ctfmodel.sqV; l--;}
      else                  {p[16]=0;            l--;}
   if (ia<=17 && l>0)
      if (astigmatic_noise) {p[17]=ctfmodel.sqrt_angle; l--;}
      else                  {p[17]=0;                   l--;}
   if (ia<=18 && l>0) {p[18]=ctfmodel.gaussian_K; l--;}
   if (ia<=19 && l>0) {p[19]=ctfmodel.sigmaU; l--;}
   if (ia<=20 && l>0)
      if (astigmatic_noise) {p[20]=ctfmodel.sigmaV; l--;}
      else                  {p[20]=0;               l--;}
   if (ia<=21 && l>0)
      if (astigmatic_noise) {p[21]=ctfmodel.gaussian_angle; l--;}
      else                  {p[21]=0;                       l--;}
   if (ia<=22 && l>0) {p[22]=ctfmodel.cU; l--;}
   if (ia<=23 && l>0)
      if (astigmatic_noise) {p[23]=ctfmodel.cV; l--;}
      else                  {p[23]=0;           l--;}
}

#define COPY_ctfmodel_TO_CURRENT_GUESS \
   assign_parameters_from_CTF(global_ctfmodel, \
      VEC_ARRAY(*global_adjust),0,CTF_PARAMETERS, \
      global_prm->astigmatic_noise);

/* Read parameters --------------------------------------------------------- */
void Adjust_CTF_Parameters::read(const FileName &fn_param) _THROW {
   FILE *fh_param; 
   if ((fh_param = fopen(fn_param.c_str(), "r")) == NULL) 
   	 REPORT_ERROR(1,(string)"Prog_Adjust_CTF::read: There is a problem " 
            "opening the file "+fn_param); 
   
   fn_ctf=get_param(fh_param,"ctf",0,"");
   fn_similar_model=get_param(fh_param,"similar_model",0,"");
   
   show_optimization=check_param(fh_param,"show_optimization"); 
   min_freq=AtoF(get_param(fh_param,"min_freq",0,"0.05")); 
   max_freq=AtoF(get_param(fh_param,"max_freq",0,"0.35")); 
   astigmatic_noise=!check_param(fh_param,"radial_noise");
   defocus_range=AtoF(get_param(fh_param,"defocus_range",0,"8000")); 
   
   particle_horizontal=AtoI(get_param(fh_param,"particle_horizontal",0,"-1"));
   particle_vertical  =AtoI(get_param(fh_param,"particle_vertical",  0,"-1"));
   if (particle_vertical==-1) particle_vertical=particle_horizontal;

   adjust.resize(CTF_PARAMETERS);
   initial_ctfmodel.enable_CTF=initial_ctfmodel.enable_CTFnoise=true;
   if (fn_similar_model=="") initial_ctfmodel.read(fn_param,FALSE);
   else                      initial_ctfmodel.read(fn_similar_model,FALSE);
   Tm=initial_ctfmodel.Tm; 
   assign_parameters_from_CTF(initial_ctfmodel,VEC_ARRAY(adjust),0,
      CTF_PARAMETERS, TRUE); 
}

/* Write to a file --------------------------------------------------------- */
void Adjust_CTF_Parameters::write(const FileName &fn_prm, bool rewrite)
   _THROW {
   ofstream fh_param;
   if (!rewrite) fh_param.open(fn_prm.c_str(),ios::app);
   else          fh_param.open(fn_prm.c_str(),ios::out);
   if (!fh_param)
      REPORT_ERROR(1,(string)"Adjust_CTF_Parameters::write: There is a problem "
            "opening the file "+fn_prm);
   fh_param << "# Adjust CTF parameters\n";
   if (fn_ctf!="")
      fh_param << "ctf="               << fn_ctf                  << endl;
   if (fn_similar_model!="")
      fh_param << "similar_model="     << fn_similar_model        << endl;
   fh_param << "min_freq="             << min_freq                << endl
            << "max_freq="             << max_freq                << endl;
   fh_param << "defocus_range="        << defocus_range           << endl;
   if (show_optimization)  fh_param    << "show_optimization=yes\n";
   if (!astigmatic_noise)  fh_param    << "radial_noise=yes\n";
   fh_param << "particle_horizontal="  << particle_horizontal     << endl;
   fh_param << "particle_vertical="    << particle_vertical       << endl;

   fh_param << initial_ctfmodel << endl;
   fh_param.close();
}

/* Show -------------------------------------------------------------------- */
void Adjust_CTF_Parameters::show() {
   cout << "CTF file:           " << fn_ctf              << endl
        << "Similar model:      " << fn_similar_model    << endl
	<< "Min Freq.:          " << min_freq            << endl
	<< "Max Freq.:          " << max_freq            << endl
	<< "Sampling:           " << Tm                  << endl
        << "Radial noise:       " << !astigmatic_noise   << endl
        << "Defocus range:      " << defocus_range       << endl
        << "Particle horizontal:" << particle_horizontal << endl
        << "Particle vertical:  " << particle_vertical   << endl
	<< "Starting:\n"          << initial_ctfmodel    << endl
    ;
}

/* Usage ------------------------------------------------------------------- */
void Adjust_CTF_Parameters::Usage() {
   cerr << "This program tries to adjust a parametric model to a CTF file.\n"
        << "Usage: adjust_ctf -i <parameters file>\n"
        << "   Where the parameters file may contain the description of a\n"
        << "      CTF and CTFnoise plus any of the following parameters\n"
	<< "   [ctf=<Fourier Xmipp Image>] : CTF file\n"
	<< "   [similar_model=<CTF and CTFnoisemodel>]: If known\n"
        << "   [show_optimization=yes]     : Show optimization process\n"
	<< "   [min_freq=<f=0.05>]         : Minimum digital frequency to use in adjust. Its value\n"
	<< "                                 should be a little lower than the dig. freq. of the first \n"  
	<< "                                 CTF zero.\n"  
 	<< "   [max_freq=<f=0.35>]         : Maximum digital frequency to use in adjust.\n"
	<< "                                 It should be higher than the last zero of the CTF.\n"  
        << "   [radial_noise=yes|no]       : By default, noise is astigmatic\n"
        << "   [particle_horizontal=-1]    : By default, the same X size as the PSD\n"
        << "   [particle_vertical=-1]      : By default, the same Y size as the PSD\n"
   ;
}

/* Produce side information ------------------------------------------------ */
void Adjust_CTF_Parameters::produce_side_info() {
   // Read the CTF file, supposed to be the uncentered squared amplitudes
   ctftomodel.read(fn_ctf);
   f=&(ctftomodel());
   
   // Resize the frequency
   global_x_digfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);
   global_y_digfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);
   global_w_digfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);
   global_x_contfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);
   global_y_contfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);
   global_w_contfreq.init_zeros(YSIZE(*f),XSIZE(*f)/2);

   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_x_digfreq) {
      XX(idx)=j; YY(idx)=i;

      // Digital frequency
      FFT_idx2digfreq(*f, idx, freq);
      global_x_digfreq(i,j)=XX(freq);
      global_y_digfreq(i,j)=YY(freq);
      global_w_digfreq(i,j)=freq.module();
      
      // Continuous frequency
      digfreq2contfreq(freq, freq, global_prm->Tm);
      global_x_contfreq(i,j)=XX(freq);
      global_y_contfreq(i,j)=YY(freq);
      global_w_contfreq(i,j)=freq.module();
   }
}

/* Generate model so far ---------------------------------------------------- */
/* The model is taken from global_adjust and global_ctfmodel is modified */
void generate_model_so_far(ImageXmipp &I, bool apply_log=FALSE) {
   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane

   assign_CTF_from_parameters(VEC_ARRAY(*global_adjust),global_ctfmodel,
      0,CTF_PARAMETERS,global_prm->astigmatic_noise);
   global_ctfmodel.Produce_Side_Info();

   I().resize(*f);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(I()) {
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(*f, idx, freq);
      digfreq2contfreq(freq, freq, global_prm->Tm);
      
      // Decide what to save
      double param;
      if (global_action<=1)
         I(i,j)=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      else if (global_action==2) {
         double E=global_ctfmodel.CTFdamping_at(XX(freq),YY(freq));
         I(i,j)=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq))+E*E;
      } else if (global_action==3 || global_action==4) {
         double ctf=global_ctfmodel.CTFpure_at(XX(freq),YY(freq));
         I(i,j)=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq))+ctf*ctf;
      } else {
         double ctf=global_ctfmodel.CTFpure_at(XX(freq),YY(freq));
         I(i,j)=ctf;
      }
      if (apply_log) I(i,j)=10*log10(I(i,j));
   }
}

/* Save intermidiate results ----------------------------------------------- */
/* First call to generate model so far and then save the image, and a couple
   of cuts along X and Y.
   
   This function returns the fitting error.*/
void save_intermidiate_results(const FileName &fn_root, bool
   generate_profiles=true) {
   
   ofstream plotX, plotY, plot_radial;
   ImageXmipp save;
   generate_model_so_far(save,false);
   
   ImageXmipp save_ctf;
   if (global_prm->particle_horizontal==-1 &&
       global_prm->particle_vertical==-1) {
      save_ctf()=global_prm->ctftomodel();

      // Substract the background from the PSD
      FOR_ALL_ELEMENTS_IN_MATRIX2D(save_ctf()) {
         double f_x=DIRECT_MAT_ELEM(global_x_contfreq,i,j);
         if (f_x<0) {
            double f_y=DIRECT_MAT_ELEM(global_y_contfreq,i,j);
            double bg=global_ctfmodel.CTFnoise_at(f_x,f_y);
            save_ctf(i,j)-=bg;
         }
      }

      // Copy the model
      FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq)
         save_ctf(i,j)=save(i,j);
   } else
      global_prm->generate_model(
         global_prm->particle_vertical,global_prm->particle_horizontal,
         save_ctf());
   save_ctf.write(fn_root+".ctfmodel");
   
   if (!generate_profiles) return;
   plotX.open((fn_root+"X.txt").c_str());
   plotY.open((fn_root+"Y.txt").c_str());
   plot_radial.open((fn_root+"_radial.txt").c_str());
   if (!plotX)
      REPORT_ERROR(1,"save_intermidiate_results::Cannot open plot file for writing\n");
   if (!plotY)
      REPORT_ERROR(1,"save_intermidiate_results::Cannot open plot file for writing\n");
   if (!plot_radial)
      REPORT_ERROR(1,"save_intermidiate_results::Cannot open plot file for writing\n");
   plotX << "# freq_dig freq_angstrom background ctf2\n";
   plotY << "# freq_dig freq_angstrom background ctf2\n";
   plot_radial << "# freq_dig freq_angstrom background ctf2\n";

   double w;
   // Generate cut along X
   for (int i=STARTINGY(save()); i<=FINISHINGY(save())/2; i++) {
      if (!global_mask(i,0)) continue;
      plotY << global_w_digfreq(i,0) << " " << global_w_contfreq(i,0)
            << " " << save(i,0) << " "
	    << (*f)(i,0) << endl;
   }

   // Generate cut along Y
   for (int j=STARTINGX(save()); j<=FINISHINGX(save())/2; j++) {
      if (!global_mask(0,j)) continue;
      plotX << global_w_digfreq(0,j) << " " << global_w_contfreq(0,j)
            << " " << save(0,j) << " "
	    << (*f)(0,j) << endl;
   }

   // Generate radial average
   matrix1D<double> radial_CTFmodel_avg(YSIZE(save())/2);
   matrix1D<double> radial_CTFampl_avg(YSIZE(save())/2);
   matrix1D<int>    radial_N(YSIZE(save())/2);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
      if (!global_mask(i,j)) continue;
      double model2=save(i,j);

      int r=FLOOR(global_w_digfreq(i,j)*(double)YSIZE(*f));
      radial_CTFmodel_avg(r)+=model2;
      radial_CTFampl_avg(r)+=(*f)(i,j);
      radial_N(r)++;
   }

   FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_CTFmodel_avg) {
      if (radial_N(i)==0) continue;
      plot_radial << global_w_digfreq(i,0) << " " << global_w_contfreq(i,0)
                  << " " << radial_CTFmodel_avg(i)/radial_N(i) << " "
	          << radial_CTFampl_avg(i)/radial_N(i) << endl;
   }

   plotX.close();
   plotY.close();
   plot_radial.close();
}

/* Generate model at a given size ------------------------------------------ */
void Adjust_CTF_Parameters::generate_model(int Ydim, int Xdim,
   matrix2D<double> &model) {
   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane

   // The right part is the scaled PSD
   model=*f;

   // Scale and remove negative values because of the interpolation
   model.self_scale_to_size_Bspline(3,Ydim,Xdim);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(model)
      if (model(i,j)<0) model(i,j)=0;
 
   // The left part is the CTF model
   assign_CTF_from_parameters(VEC_ARRAY(*global_adjust),global_ctfmodel,
      0,CTF_PARAMETERS,global_prm->astigmatic_noise);
   global_ctfmodel.Produce_Side_Info();

   FOR_ALL_ELEMENTS_IN_MATRIX2D(model) {
      if (j>Xdim/2) continue;
      
      XX(idx)=j; YY(idx)=i;
      FFT_idx2digfreq(model, idx, freq);
      digfreq2contfreq(freq, freq, global_prm->Tm);
      
      model(i,j)=global_ctfmodel.CTFpure_at(XX(freq),YY(freq));
   }
}

/* Build mask -------------------------------------------------------------- */
void build_freq_mask() {
   global_mask.init_zeros(global_w_digfreq);
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
      if (global_value_th!=-1 && (*f)(i,j)>=global_value_th)
         continue;
      if (global_w_digfreq(i,j)>=global_max_freq ||
          global_w_digfreq(i,j)<=global_min_freq) continue;
      global_mask(i,j)=1;
   }
}

/* CTF fitness ------------------------------------------------------------- */
/* This function measures the distance between the estimated CTF and the
   measured CTF */
double CTF_fitness(double *p) {
   double retval;

   // Generate CTF model
   switch (global_action) {
      // Remind that p is a vector whose first element is at index 1
      case 0: assign_CTF_from_parameters(p-FIRST_SQRT_PARAMETER+1,
                global_ctfmodel,FIRST_SQRT_PARAMETER,SQRT_CTF_PARAMETERS,
                global_prm->astigmatic_noise);
              if (global_show>=2) {
                 cout << "Input vector:";
                 for (int i=1; i<=SQRT_CTF_PARAMETERS; i++) cout << p[i] << " ";
                 cout << endl;
              }   
              break;
      case 1: assign_CTF_from_parameters(p-FIRST_SQRT_PARAMETER+1,
                global_ctfmodel,FIRST_SQRT_PARAMETER,BACKGROUND_CTF_PARAMETERS,
                global_prm->astigmatic_noise);
              if (global_show>=2) {
                 cout << "Input vector:";
                 for (int i=1; i<=BACKGROUND_CTF_PARAMETERS; i++) cout << p[i] << " ";
                 cout << endl;
              }   
              break;
      case 2: assign_CTF_from_parameters(p-FIRST_ENVELOPE_PARAMETER+1,
                global_ctfmodel,FIRST_ENVELOPE_PARAMETER,ENVELOPE_PARAMETERS,
                global_prm->astigmatic_noise);
              if (global_show>=2) {
                 cout << "Input vector:";
                 for (int i=1; i<=ENVELOPE_PARAMETERS; i++) cout << p[i] << " ";
                 cout << endl;
              }   
              break;
      case 3: assign_CTF_from_parameters(p-FIRST_DEFOCUS_PARAMETER+1,
                global_ctfmodel,FIRST_DEFOCUS_PARAMETER,DEFOCUS_PARAMETERS,
                global_prm->astigmatic_noise);
              if (global_show>=2) {
                 cout << "Input vector:";
                 for (int i=1; i<=DEFOCUS_PARAMETERS; i++) cout << p[i] << " ";
                 cout << endl;
              }   
              break;
      case 4: assign_CTF_from_parameters(p-0+1,
                global_ctfmodel,0,CTF_PARAMETERS,
                global_prm->astigmatic_noise);
              if (global_show>=2) {
                 cout << "Input vector:";
                 for (int i=1; i<=CTF_PARAMETERS; i++) cout << p[i] << " ";
                 cout << endl;
              }   
              break;
   }
   global_ctfmodel.Produce_Side_Info();
   if (global_show>=2) cout << "Model:\n" << global_ctfmodel << endl;
   if(!global_ctfmodel.physical_meaning()) {
      if (global_show>=2) cout << "Does not have physical meaning\n";
      return global_heavy_penalization;
   }
      
   // Now the 2D error
   double distsum=0;
   int N=0;
   for (int i=0; i<YSIZE(global_w_digfreq); i+=global_evaluation_reduction)
      for (int j=0; j<XSIZE(global_w_digfreq); j+=global_evaluation_reduction) {
         if (!DIRECT_MAT_ELEM(global_mask,i,j)) continue;

         // Compute each component
         double f_x=DIRECT_MAT_ELEM(global_x_contfreq,i,j);
         double f_y=DIRECT_MAT_ELEM(global_y_contfreq,i,j);
         double bg=global_ctfmodel.CTFnoise_at(f_x,f_y);
         double envelope, ctf_without_damping, ctf_with_damping;
         double ctf2_th;
         switch (global_action) {
            case 0:
            case 1: ctf2_th=bg; break;
            case 2: envelope=global_ctfmodel.CTFdamping_at(f_x,f_y);
                    ctf2_th=bg+envelope*envelope;
                    break;
            case 3:
            case 4: 
                    if (global_prm->initial_ctfmodel.DeltafU!=0) {
                       // If there is an initial model, the true solution
                       // cannot be too far
                       if (ABS(global_prm->initial_ctfmodel.DeltafU-
                               global_ctfmodel.DeltafU)>10000 ||
                           ABS(global_prm->initial_ctfmodel.DeltafV-
                               global_ctfmodel.DeltafV)>10000) {
                          if (global_show>=2)
                             cout << "Too far from hint\n";
                          return global_heavy_penalization;
                       }
                    }
                    envelope=global_ctfmodel.CTFdamping_at(f_x,f_y);
                    ctf_without_damping=global_ctfmodel.CTFpure_without_damping_at(f_x,f_y);
                    ctf_with_damping=envelope*ctf_without_damping;
                    ctf2_th=bg+ctf_with_damping*ctf_with_damping;
                    break;
         }
         
         // Compute distance
         double ctf2=DIRECT_MAT_ELEM(*f,i,j);
         double dist;
         switch (global_action) {
            case 0:
            case 1:
               dist=ABS(ctf2-bg);
               if (global_penalize && bg>ctf2 &&
                   DIRECT_MAT_ELEM(global_w_digfreq,i,j)>global_max_gauss_freq)
                  dist*=global_current_penalty;
               break;
            case 2:
               dist=ABS(ctf2-ctf2_th);
               if (global_penalize && ctf2_th<ctf2 &&
                   DIRECT_MAT_ELEM(global_w_digfreq,i,j)>global_max_gauss_freq)
                  dist*=global_current_penalty;
               break;
            case 3:
            case 4:
               if (envelope>XMIPP_EQUAL_ACCURACY)
                  dist=ABS(ctf2-ctf2_th)/(envelope*envelope);
               else dist=ABS(ctf2-ctf2_th);
                  // This expression comes from mapping any value so that
                  // bg becomes 0, and bg+envelope^2 becomes 1
                  // This is the transformation
                  //        (x-bg)      x-bg
                  //    -------------=-------
                  //    (bg+env^2-bg)  env^2
                  // If we substract two of this scaled values
                  //    x-bg      y-bg       x-y
                  //   ------- - ------- = ------- 
                  //    env^2     env^2     env^2
               break;
         }
         distsum+=dist;
         N++;
      }
   if (N>0) retval=distsum/N;
   else     retval=global_heavy_penalization;

   // Show some debugging information
   if (global_show>=2) {
      cout << "Fitness=" << retval << endl;
      if (global_show==3) {
         save_intermidiate_results("PPP");
         cout << "Press any key\n"; char c; cin >> c;
      }
   }

   return retval;
}

/* Compute central region -------------------------------------------------- */
void compute_central_region(double &w1, double &w2, double ang) {
   w1=global_max_freq; w2=global_min_freq; 
   matrix1D<double> freq(2), dir(2);
   
   // Compute first and third zero in the given direction
   VECTOR_R2(dir,COSD(ang),SIND(ang));

   // Detect first zero
   global_ctfmodel.zero(1,dir,freq);
   if (XX(freq)==-1 && YY(freq)==-1) w1=global_min_freq;
   else {
      contfreq2digfreq(freq,freq,global_prm->Tm);
      double w;
      if (XX(dir)>0.1) w=XX(freq)/XX(dir);
      else w=YY(freq)/YY(dir);
      w1=MAX(global_min_freq,MIN(w1,w));
   }

   // Detect fifth zero
   global_ctfmodel.zero(5,dir,freq);
   if (XX(freq)==-1 && YY(freq)==-1) w2=global_max_freq;
   else {
      double w;
      contfreq2digfreq(freq,freq,global_prm->Tm);
      if (XX(dir)>0.1) w=XX(freq)/XX(dir);
      else w=YY(freq)/YY(dir);
      w2=MIN(global_max_freq,MAX(w2,w));
   }
}

/* Center focus ----------------------------------------------------------- */
void center_optimization_focus(bool adjust_freq, bool adjust_th,
   double margin=1) {
   if (global_prm->show_optimization)
      cout << "Freq frame before focusing=" << global_min_freq << ","
           << global_max_freq << endl
	   << "Value_th before focusing=" << global_value_th << endl;

   double w1=global_min_freq,w2=global_max_freq;
   if (adjust_freq) {
      double w1U, w2U, w1V, w2V;
      compute_central_region(w1U,w2U,global_ctfmodel.azimuthal_angle);
      compute_central_region(w1V,w2V,global_ctfmodel.azimuthal_angle+90);
      w1=MIN(w1U,w1V);
      w2=MAX(w2U,w2V);
      global_min_freq=MAX(global_min_freq,w1-0.05);
      global_max_freq=MIN(global_max_freq,w2+0.01);
   }

   // Compute maximum value within central region
   if (adjust_th) {
      ImageXmipp save; generate_model_so_far(save);
      double max_val=0;
      FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
         double w=global_w_digfreq(i,j);
         if (w>=w1 && w<=w2)
   	    max_val=MAX(max_val,save(i,j));
      }      
      if (global_value_th!=-1)
         global_value_th=MIN(global_value_th,max_val*margin);
      else global_value_th=max_val*margin;
   }

   if (global_prm->show_optimization)
      cout << "Freq frame after focusing=" << global_min_freq << ","
           << global_max_freq << endl
	   << "Value_th after focusing=" << global_value_th << endl;
}

// Estimate sqrt parameters ------------------------------------------------
// Results are written in global_ctfmodel
void estimate_background_sqrt_parameters() {
   if (global_prm->show_optimization)
      cout << "Computing first sqrt background ...\n";

   // Estimate the base line taking the value of the CTF
   // for the maximum X and Y frequencies
   double base_line=0; int N=0;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq)
      if (global_w_digfreq(i,j)>0.4) {N++; base_line+=(*f)(i,j);}
   global_ctfmodel.base_line=base_line/N;
   
   // Find the linear least squares solution for the sqrt part
   matrix2D<double> A(2,2); A.init_zeros();
   matrix1D<double> b(2);   b.init_zeros();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
      if (!global_mask(i,j)) continue;
      
      // Compute weight for this point
      double weight=1+global_max_freq-global_w_digfreq(i,j);

      // Compute error
      double explained=global_ctfmodel.CTFnoise_at(
         global_x_contfreq(i,j),global_y_contfreq(i,j));
      double unexplained=(*f)(i,j)-explained;
      if (unexplained<=0) continue;
      unexplained=log(unexplained);

      double X=-sqrt(global_w_contfreq(i,j));
      A(0,0)+=weight*X*X;
      A(0,1)+=weight*X;
      A(1,1)+=weight*1;
      b(0)  +=X*weight*unexplained;
      b(1)  +=  weight*unexplained;
   }
   A(1,0)=A(0,1);
   b=A.inv()*b;
   global_ctfmodel.sqU=global_ctfmodel.sqV=b(0);
   global_ctfmodel.sqrt_K=exp(b(1));
   global_ctfmodel.sqrt_angle=0;
   
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "First SQRT Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step01_first_sqrt_fit");
   }
   
   // Now optimize .........................................................
   double fitness;
   matrix1D<double> steps;
   steps.resize(SQRT_CTF_PARAMETERS);
   steps.init_constant(0);
   steps(0)=steps(1)=steps(2)=1;
   if (global_prm->astigmatic_noise) {steps(3)=steps(4)=1;}

   // Optimize without penalization
   if (global_prm->show_optimization)
      cout << "Looking for best fitting sqrt ...\n";
   global_penalize=false;
   int iter;
   Powell_optimizer(*global_adjust, FIRST_SQRT_PARAMETER+1, SQRT_CTF_PARAMETERS,
      &CTF_fitness, 0.05, fitness, iter, steps,
      global_prm->show_optimization);

   // Optimize with penalization
   if (global_prm->show_optimization)
      cout << "Penalizing best fitting sqrt ...\n";
   global_penalize=true;
   global_current_penalty=2;
   int imax=CEIL(log(global_penalty)/log(2.0));
   for (int i=1; i<=imax; i++) {
      if (global_prm->show_optimization)
         cout << "     Iteration " << i
              << " penalty=" << global_current_penalty << endl;
      Powell_optimizer(*global_adjust, FIRST_SQRT_PARAMETER+1,
         SQRT_CTF_PARAMETERS, &CTF_fitness,
         0.05, fitness, iter, steps, global_prm->show_optimization);
      global_current_penalty*=2;
      global_current_penalty=MIN(global_current_penalty, global_penalty);
   }
   // Keep the result in global_prm->adjust
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "Best penalized SQRT Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step02_best_penalized_sqrt_fit");
   }

   center_optimization_focus(false,true,1.5);
}

// Estimate gaussian parameters --------------------------------------------
//#define DEBUG
void estimate_background_gauss_parameters() {
   if (global_prm->show_optimization)
      cout << "Computing first background parameters ...\n";

   // Compute radial averages
   matrix1D<double> radial_CTFmodel_avg(YSIZE(*f)/2);
   matrix1D<double> radial_CTFampl_avg(YSIZE(*f)/2);
   matrix1D<int>    radial_N(YSIZE(*f)/2);
   double w_max_gauss=0.25;
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
      if (!global_mask(i,j)) continue;
      double w=global_w_digfreq(i,j);
      if (w>w_max_gauss) continue;

      int r=FLOOR(w*(double)YSIZE(*f));
      radial_CTFmodel_avg(r)+=global_ctfmodel.CTFnoise_at(
         global_x_contfreq(i,j),global_y_contfreq(i,j));
      radial_CTFampl_avg(r)+=(*f)(i,j);
      radial_N(r)++;
   }

   // Compute the average radial error
   double error2_avg=0;
   int N_avg=0;
   matrix1D<double> error; error.init_zeros(radial_CTFmodel_avg);
   FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_CTFmodel_avg) {
      if (radial_N(i)==0) continue;
      error(i)=(radial_CTFampl_avg(i)-radial_CTFmodel_avg(i))/radial_N(i);
      error2_avg+=error(i)*error(i); N_avg++;
   }
   if (N_avg!=0) error2_avg/=N_avg;
   
   #ifdef DEBUG
      cout << "Error2 avg=" << error2_avg << endl;
   #endif

   // Compute the minimum radial error
   bool   first=true, OK_to_proceed=false;
   double error2_min=0, wmin, fmin;
   FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_CTFmodel_avg) {
      if (radial_N(i)==0) continue;
      double w=global_w_digfreq(i,0);

      if (error(i)<0 && first) continue;
      else if (error(i)<0) break;
      double error2=error(i)*error(i);
      // If the two lines cross, do not consider any error until
      // the cross is "old" enough
      if (first && error2>0.15*error2_avg) OK_to_proceed=true;
      if (first && i>0) OK_to_proceed&=(error(i)<error(i-1));

      // If the error now is bigger than a 30% (1.69=1.3*1.3) of the error min
      // this must be a rebound. Stop here
      if (!first && error2>1.69*error2_min) break;
      if ( first && OK_to_proceed) {wmin=w; error2_min=error2; first=false;}
      if (!first && error2<error2_min) {wmin=w; error2_min=error2;}
      #ifdef DEBUG
	 cout << w << " " << error2 << " " << wmin << " " << endl;
      #endif
   }
   
   // Compute the frequency of the minimum error
   global_max_gauss_freq=wmin;
   fmin=wmin/global_prm->Tm;
   #ifdef DEBUG
      cout << "Freq of the minimum error: " << wmin << " " << fmin << endl;
   #endif

   // Compute the maximum radial error
   first=true;
   double error2_max=0, wmax, fmax;
   FOR_ALL_ELEMENTS_IN_MATRIX1D(radial_CTFmodel_avg) {
      if (radial_N(i)==0) continue;
      double w=global_w_digfreq(i,0);
      if (w>wmin) continue;

      if (error(i)<0 && first) continue;
      else if (error(i)<0) break;
      double error2=error(i)*error(i);
      if (first) {wmax=w; error2_max=error2; first=false;}
      if (error2>error2_max) {wmax=w; error2_max=error2;}
      #ifdef DEBUG
	 cout << w << " " << error2 << " " << wmax << endl;
      #endif
   }
   fmax=global_ctfmodel.cV=global_ctfmodel.cU=wmax/global_prm->Tm;
   #ifdef DEBUG
      cout << "Freq of the maximum error: " << wmax << " " << fmax << endl;
   #endif

   // Find the linear least squares solution for the gauss part
   matrix2D<double> A(2,2); A.init_zeros();
   matrix1D<double> b(2);   b.init_zeros();
   FOR_ALL_ELEMENTS_IN_MATRIX2D(global_w_digfreq) {
      if (!global_mask(i,j)) continue;
      if (global_w_digfreq(i,j)>wmin) continue;
      double fmod=global_w_contfreq(i,j);

      // Compute weight for this point
      double weight=1+global_max_freq-global_w_digfreq(i,j);

      // Compute error
      double explained=global_ctfmodel.CTFnoise_at(
         global_x_contfreq(i,j),global_y_contfreq(i,j));
      double unexplained=(*f)(i,j)-explained;
      if (unexplained<=0) continue;
      unexplained=log(unexplained);
      double F=-(fmod-fmax)*(fmod-fmax);
      A(0,0)+=weight*1;
      A(0,1)+=weight*F;
      A(1,1)+=weight*F*F;
      b(0)  +=  weight*unexplained;
      b(1)  +=F*weight*unexplained;
   }
   A(1,0)=A(0,1);
   b=A.inv()*b;
   global_ctfmodel.sigmaU=MIN(ABS(b(1)),95e3); // This value should be
   global_ctfmodel.sigmaV=MIN(ABS(b(1)),95e3); // conformant with the physical
      	             	      	               // meaning routine in CTF.cc
   global_ctfmodel.gaussian_K=exp(b(0));

   // Store the CTF values in global_prm->adjust
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "First Background Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step03_first_background_fit");
   }

   center_optimization_focus(false,true,1.5);
}
#undef DEBUG

// Estimate envelope parameters --------------------------------------------
//#define DEBUG
void estimate_envelope_parameters() {
   if (global_prm->show_optimization)
      cout << "Looking for best fitting envelope ...\n";

   // Set the envelope
   global_ctfmodel.Ca              = 2.0; // Typical value
   global_ctfmodel.K               = 1.0;
   global_ctfmodel.espr            = 0.0;
   global_ctfmodel.ispr            = 0.0;
   global_ctfmodel.alpha           = 0.0;
   global_ctfmodel.DeltaF          = 0.0;
   global_ctfmodel.DeltaR          = 0.0;
   global_ctfmodel.Q0              = 0.0;
   COPY_ctfmodel_TO_CURRENT_GUESS;

   // Now optimize the envelope
   global_penalize=false;
   int iter;
   double fitness;
   matrix1D<double> steps;
   steps.resize(ENVELOPE_PARAMETERS);
   steps.init_constant(1);
   steps(1)=0; // Do not optimize Cs
   steps(5)=0; // Do not optimize for alpha, since Ealpha depends on the
               // defocus
   Powell_optimizer(*global_adjust, FIRST_ENVELOPE_PARAMETER+1, ENVELOPE_PARAMETERS,
      &CTF_fitness, 0.05, fitness, iter, steps,
      global_prm->show_optimization);

   // Keep the result in global_prm->adjust
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "Best envelope Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step05_best_envelope_fit");
   }

   // Optimize with penalization
   if (global_prm->show_optimization)
      cout << "Penalizing best fitting envelope ...\n";
   global_penalize=true;
   global_current_penalty=2;
   int imax=CEIL(log(global_penalty)/log(2.0));
   for (int i=1; i<=imax; i++) {
      if (global_prm->show_optimization)
         cout << "     Iteration " << i
              << " penalty=" << global_current_penalty << endl;
      Powell_optimizer(*global_adjust, FIRST_ENVELOPE_PARAMETER+1,
         ENVELOPE_PARAMETERS, &CTF_fitness,
         0.05, fitness, iter, steps, global_prm->show_optimization);
      global_current_penalty*=2;
      global_current_penalty=MIN(global_current_penalty, global_penalty);
   }
   // Keep the result in global_prm->adjust
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "Best envelope Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step06_best_penalized_envelope_fit");
   }
}
#undef DEBUG

// Estimate  defoci --------------------------------------------------------
//#define DEBUG
void estimate_defoci() {
   if (global_prm->show_optimization)
      cout << "Looking for first defoci ...\n";
   double best_defocusU, best_defocusV, best_error, best_angle;
   bool first=TRUE;
   int i,j;
   double defocusV, defocusU;

   double defocusV0=-1e3, defocusU0=-1e3;
   double defocusVF=-100e3, defocusUF=-100e3;
   double initial_defocusStep=16e3;
   matrix2D<double> error;
   
   // Check if there is no initial guess
   double min_allowed_defocusU=-100e3, max_allowed_defocusU=-1e3;
   double min_allowed_defocusV=-100e3, max_allowed_defocusV=-1e3;
   if (global_prm->initial_ctfmodel.DeltafU!=0) {
      initial_defocusStep=global_prm->defocus_range;
      max_allowed_defocusU=defocusU0=
         MIN(-1000,global_prm->initial_ctfmodel.DeltafU+global_prm->defocus_range);
      min_allowed_defocusU=defocusUF=
         MAX(-100000,global_prm->initial_ctfmodel.DeltafU-global_prm->defocus_range);
      if (global_prm->initial_ctfmodel.DeltafV==0) {
         defocusV0=defocusU0; min_allowed_defocusV=min_allowed_defocusU;
         defocusVF=defocusUF; max_allowed_defocusV=max_allowed_defocusU;
      } else {
         max_allowed_defocusV=defocusV0=
            MIN(-1000,global_prm->initial_ctfmodel.DeltafV+global_prm->defocus_range);
         min_allowed_defocusV=defocusVF=
            MAX(-100000,global_prm->initial_ctfmodel.DeltafV-global_prm->defocus_range);
      }
   }
   
   matrix1D<double> steps(DEFOCUS_PARAMETERS);
   steps.init_constant(1);
   for (double defocusStep=initial_defocusStep; defocusStep>=global_prm->defocus_range; defocusStep/=2) {
      error.resize(CEIL((defocusV0-defocusVF)/defocusStep+1),
                   CEIL((defocusU0-defocusUF)/defocusStep+1));
      if (global_prm->show_optimization)
         cout << "V=["<<defocusV0 << "," << defocusVF << "]\n"
              << "U=["<<defocusU0 << "," << defocusUF << "]\n"
              << "Defocus step=" << defocusStep << endl;
      for (defocusV=defocusV0,i=0; defocusV>=defocusVF; defocusV-=defocusStep, i++) {
         for (defocusU=defocusU0,j=0; defocusU>=defocusUF; defocusU-=defocusStep, j++) {
            bool first_angle=true;
            for (double angle=0; angle<90; angle+=45) {
               int iter;
               double fitness;

               (*global_adjust)(0)=defocusU;
               (*global_adjust)(1)=defocusV;
               (*global_adjust)(2)=angle;
            
               Powell_optimizer(*global_adjust, FIRST_DEFOCUS_PARAMETER+1,
                  DEFOCUS_PARAMETERS, &CTF_fitness,
                  0.05, fitness, iter, steps, false);
               
               if ((first_angle || fitness<error(i,j)) &&
                   (global_ctfmodel.DeltafU>=min_allowed_defocusU &&
                    global_ctfmodel.DeltafU<=max_allowed_defocusU &&
                    global_ctfmodel.DeltafV>=min_allowed_defocusV &&
                    global_ctfmodel.DeltafV<=max_allowed_defocusV)) {
                  error(i,j)=fitness;
                  first_angle=false;
                  if (error(i,j)<best_error || first) {
                     best_error=error(i,j);
                     best_defocusU=global_ctfmodel.DeltafU;
                     best_defocusV=global_ctfmodel.DeltafV;
                     best_angle=global_ctfmodel.azimuthal_angle;
                     first=false;
                     if (global_prm->show_optimization)
                        cout << "    (DefocusU,DefocusV)=(" << defocusU << ","
                             << defocusV << "), ang=" << angle
                             << " --> (" << global_ctfmodel.DeltafU << ","
                             << global_ctfmodel.DeltafV << "),"
                             << global_ctfmodel.azimuthal_angle << " error="
                             << error(i,j) << endl;
                  }
               }
            }
         }
      }

      // Compute the range of the errors
      double errmin=error(0,0), errmax=error(0,0);
      for (int ii=STARTINGY(error); ii<=FINISHINGY(error); ii++)
         for (int jj=STARTINGX(error); jj<=FINISHINGX(error); jj++) {
            if (error(ii,jj)!=global_heavy_penalization)
               if (error(ii,jj)<errmin) errmin=error(ii,jj);
               else if (errmax==global_heavy_penalization) errmax=error(ii,jj);
               else if (error(ii,jj)>errmax) errmax=error(ii,jj);
         }

      // Find those defoci which are within a 10% of the maximum
      if (global_show>=2) cout << "Range=" << errmax-errmin << endl;
      double best_defocusVmin=best_defocusV, best_defocusVmax=best_defocusV;
      double best_defocusUmin=best_defocusU, best_defocusUmax=best_defocusU;
      for (defocusV=defocusV0,i=0; defocusV>defocusVF; defocusV-=defocusStep, i++) {
         for (defocusU=defocusU0,j=0; defocusU>defocusUF; defocusU-=defocusStep, j++) {
           if (global_show>=2)
               cout << i << "," << j << " " << error(i,j) << " " << defocusU << " " << defocusV << endl
                    << best_defocusUmin << " " << best_defocusUmax << endl
                    << best_defocusVmin << " " << best_defocusVmax << endl;
            if ((error(i,j)-errmin)/(errmax-errmin)<=0.1) {
               if (defocusV<best_defocusVmin) best_defocusVmin=defocusV;
               if (defocusU<best_defocusUmin) best_defocusUmin=defocusU;
               if (defocusV>best_defocusVmax) best_defocusVmax=defocusV;
               if (defocusU>best_defocusUmax) best_defocusUmax=defocusU;
            }
         }
      }

      defocusVF=MAX(min_allowed_defocusV,best_defocusVmin-defocusStep);
      defocusV0=MIN(max_allowed_defocusV,best_defocusVmax+defocusStep);
      defocusUF=MAX(min_allowed_defocusU,best_defocusUmin-defocusStep);
      defocusU0=MIN(max_allowed_defocusU,best_defocusUmax+defocusStep);
      i=j=0;
      if (global_show>=2) {
         ImageXmipp save; save()=error; save.write("error.xmp");
         cout << "Press any key: Error saved\n";
         char c; cin >> c;
      }
   }

   global_ctfmodel.DeltafU=best_defocusU;
   global_ctfmodel.DeltafV=best_defocusV;
   global_ctfmodel.azimuthal_angle=best_angle;

   // Keep the result in global_prm->adjust
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "First defocus Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step07_first_defocus_fit");
   }
}
#undef DEBUG

void test() {
   cout << "Test\n";
   (*global_adjust)(0)=-1000;
   (*global_adjust)(1)=-1000;
   (*global_adjust)(2)=0;

   assign_CTF_from_parameters(VEC_ARRAY(*global_adjust),global_ctfmodel,
      0,CTF_PARAMETERS,global_prm->astigmatic_noise);
   global_ctfmodel.Produce_Side_Info();

   cout << global_ctfmodel << endl;

   matrix1D<int>    idx(2);  // Indexes for Fourier plane
   matrix1D<double> freq(2); // Frequencies for Fourier plane
   for (int j=0; j<XSIZE(*f)/2; j++) {
      if (!global_mask(0,j)) continue;
      XX(idx)=j; YY(idx)=0;
      FFT_idx2digfreq(*f, idx, freq);
      digfreq2contfreq(freq, freq, global_prm->Tm);
      
      double bg=global_ctfmodel.CTFnoise_at(XX(freq),YY(freq));
      double E=global_ctfmodel.CTFdamping_at(XX(freq),YY(freq),true);
      double ctf=global_ctfmodel.CTFpure_without_damping_at(XX(freq),YY(freq),true);
      cout << global_w_digfreq(0,j) << " "
           << bg << " " << bg+E*E << " " << ctf << " "
           << (*f)(0,j) << endl;
   }

   save_intermidiate_results("PPP");
   cout << "Press any key\n"; char c; cin >> c;
}

/* Main routine ------------------------------------------------------------ */
//#define DEBUG
double ROUT_Adjust_CTF(Adjust_CTF_Parameters &prm, bool standalone) {
   global_prm=&prm;
   if (standalone || prm.show_optimization) prm.show();
   prm.produce_side_info();

   // Set some global variables
   global_adjust=&prm.adjust;
   global_penalize=false;
   global_max_gauss_freq=0;
   global_heavy_penalization=f->compute_max()*XSIZE(*f)*YSIZE(*f);
   global_show=0;

   // Build initial frequency mask
   global_value_th=-1;
   global_min_freq=prm.min_freq;
   global_max_freq=prm.max_freq;
   build_freq_mask();

   // Some variables needed by all steps
   int iter;
   double fitness;
   matrix1D<double> steps;

   /************************************************************************
     STEPs 0 & 1:  Find background which best fits the CTF
   /************************************************************************/
   global_ctfmodel.enable_CTFnoise=true;
   global_ctfmodel.enable_CTF=false;
   global_evaluation_reduction=4;
   
   // If initial parameters weren�t supplied for the gaussian curve,
   // estimate them from the CTF file
   global_action=0;
   if (prm.adjust(FIRST_SQRT_PARAMETER)==0) {
      estimate_background_sqrt_parameters();
      estimate_background_gauss_parameters();
   }

   // Optimize the current background
   global_action=1;
   global_penalize=true;
   global_current_penalty=global_penalty;
   steps.resize(BACKGROUND_CTF_PARAMETERS);
   steps.init_constant(1);
   if (!global_prm->astigmatic_noise)
      steps(3)=steps(4)=steps(7)=steps(8)=steps(10)=0;
   Powell_optimizer(*global_adjust, FIRST_SQRT_PARAMETER+1,
      BACKGROUND_CTF_PARAMETERS, &CTF_fitness,
      0.01, fitness, iter, steps, global_prm->show_optimization);
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "Best background Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step04_best_background_fit");
   }

   /************************************************************************
     STEP 2:  Find envelope which best fits the CTF
   /************************************************************************/
   global_action=2;
   global_ctfmodel.enable_CTF=true;
   if (prm.initial_ctfmodel.K==0) {
      global_ctfmodel.kV              = prm.initial_ctfmodel.kV;
      global_ctfmodel.Cs              = prm.initial_ctfmodel.Cs;
      estimate_envelope_parameters();
   } else {
      global_ctfmodel.K               = prm.initial_ctfmodel.K;
      global_ctfmodel.kV              = prm.initial_ctfmodel.kV;
      global_ctfmodel.DeltafU         = prm.initial_ctfmodel.DeltafU;
      global_ctfmodel.DeltafV         = prm.initial_ctfmodel.DeltafV;
      global_ctfmodel.azimuthal_angle = prm.initial_ctfmodel.azimuthal_angle;
      global_ctfmodel.Cs              = prm.initial_ctfmodel.Cs;
      global_ctfmodel.Ca              = prm.initial_ctfmodel.Ca;
      global_ctfmodel.espr            = prm.initial_ctfmodel.espr;
      global_ctfmodel.ispr            = prm.initial_ctfmodel.ispr;
      global_ctfmodel.alpha           = prm.initial_ctfmodel.alpha;
      global_ctfmodel.DeltaF          = prm.initial_ctfmodel.DeltaF;
      global_ctfmodel.DeltaR          = prm.initial_ctfmodel.DeltaR;
      global_ctfmodel.Q0              = prm.initial_ctfmodel.Q0;
      COPY_ctfmodel_TO_CURRENT_GUESS;
   }

   /************************************************************************
     STEP 3:  the defocus and angular parameters
   /************************************************************************/
   global_action=3;
   estimate_defoci();

   /************************************************************************
     STEP 4:  all parameters
   /************************************************************************/
   global_action=4;

   steps.resize(CTF_PARAMETERS);
   steps.init_constant(1);
   steps(3)=0; // kV
   steps(5)=0; // The spherical aberration (Cs) is not optimized
   if (!global_prm->astigmatic_noise)
      steps(16)=steps(17)=steps(20)=steps(21)=steps(23)=0;
   Powell_optimizer(*global_adjust, 0+1,
      CTF_PARAMETERS, &CTF_fitness,
      0.01, fitness, iter, steps, global_prm->show_optimization);
   COPY_ctfmodel_TO_CURRENT_GUESS;

   if (global_prm->show_optimization) {
      cout << "Best fast Fit:\n" << global_ctfmodel << endl;
      save_intermidiate_results("step08_best_fast_fit");
   }
   
   #ifdef NEVER_DEFINED
      global_evaluation_reduction=1;
      Powell_optimizer(*global_adjust, 0+1,
         CTF_PARAMETERS, &CTF_fitness,
         0.01, fitness, iter, steps, global_prm->show_optimization);
      COPY_ctfmodel_TO_CURRENT_GUESS;

      if (global_prm->show_optimization) {
         cout << "Best fast Fit:\n" << global_ctfmodel << endl;
         save_intermidiate_results("step09_best_fit");
      }
   #endif
   
   /************************************************************************
     STEP 5:  Produce output
   /************************************************************************/
   FileName fn_root=prm.fn_ctf.without_extension();
   global_action=5;
   save_intermidiate_results(fn_root,false);
   global_ctfmodel.write(fn_root+".ctfparam");
   
   return 0.0;
}
