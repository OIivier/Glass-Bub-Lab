
/*
 * FHN model simulation with diffusion using the finite differences method.
 * 
 * UPDATED:
 * - Transforms the 2D spatial loop in 1D loop by using the positions 
 *	pointed by *PosX and *PosY. 
 * 
 * - Exports to matlab.
 */

#include "mex.h"
#include "omp.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // exp
#include <cstdlib>
#include <ctime>

using namespace std;

// Definitions
#define SKIP_VDATA 1		// skip points for video 
#define SKIP_VTS   1		// skip points for time series


// Structures and classes ------------------------------------------- // 
// Struct to count the number of objects created

// Class for pacemakers
class Pacemaker 
{
  public:
    double wLP;
    double Ip;
    int x;
    int y;
    int r;
    int on;
    int off;
};


//int main(){ 

void myspiral(mwSize Nt, double dx, double dt, double D, double wH, double wLP, 
		    mwSize pon, mwSize poff, mwSize pdx, mwSize pdy, mwSize TotalBeats, 
		    mwSize nPosXY, 
		    double *PosX, double *PosY, int NoiseSignal, int NoiseRecord, mwSize nPosXYMod,
		    double *PosXMod, double *PosYMod, mwSize MapModTime,
		    double *outVid, double *outPace)  
{
  
  // Grid parameters
  mwSize Nx = 82; //Width of grid in pixels
  mwSize Ny = Nx; //Height of grid in pixels, same as width
  mwSize xHalf = Nx/2; //Half of grid width
  mwSize yHalf = Ny/2; //Half of grid height

  // Diffusion and integration steps
//   double dx = 0.025;
//   double dt = 0.025;	// 2-dimentional: (explicit scheme) convergent if 2*D*dt/(dx^2) <= 1/2;
  
  // FHN parameters - FitzHugh–Nagumo equations, see paper; these are constants
  double wh    = wH; //Obtained from previous paper, value set in MATLAB file
  double wL    = 0.4;//0.4; //Obtained from previous paper
  double beta  = 0.7;//0.7; //Obtained from previous paper
  double gamma = 0.5;//0.5; //Obtained from previous paper
  double eps   = 0.44;//0.44; //0.44//Obtained from previous paper
  
  // Monolayer radius
//   mwSize LayerRadius = (Nx-2)/2;
  
  // Threshold for spikes
  double threshold = 0.25;
  int CountBeats = 0; //Beat counter
  
  // Channel obstacle parameters ------------------------------------ //
  // Ellipse
  mwSize Cx  = xHalf, Cy  = yHalf;	// center
  // ---------------------------------------------------------------- //
  
  
  // Pacemaker and Probe parameters --------------------------------- //
  vector<Pacemaker> P;
  
  P.resize(1);
  
  P[0].wLP = wLP;			// affects frequency of pacing 
  P[0].Ip = 2;				// pacing current 
  P[0].x = Cx + pdx;			// x-coordinate
  P[0].y = Cy + pdy;			// y-coordinate
  P[0].r = 4; 	 			// pacemaker radius
  P[0].on = pon;			// time at which this pacemaker starts pacing
  P[0].off= poff;			// time at which this pacemaker stops  pacing
    
//   P[1].wLP = 0.1;			// affects frequency of pacing
//   P[1].Ip = 2;				// pacing current
//   P[1].x = 68;        // x-coordinate
//   P[1].y = 48;	    // y-coordinate
//   P[1].r = 4; 	 			// pacemaker radius
//   P[1].on = 14000;			// time at which this pacemaker starts pacing
//   P[1].off= 20000;			// time at which this pacemaker stops  pacing
//   


 
  // ---------------------------------------------------------------- //
  
  vector<int> posx(nPosXY);
  vector<int> posy(nPosXY);
  
   #pragma omp parallel for //Parallelize what's inside
    {
      for(int i = 0; i < nPosXY; ++i){
	    posx[i] = (int) PosX[i];
	    posy[i] = (int) PosY[i];  
      }
    }
    
  vector<int> posxMod(nPosXYMod);
  vector<int> posyMod(nPosXYMod);
  
   #pragma omp parallel for //Parallelize what's inside
    {
      for(int i = 0; i < nPosXYMod; ++i){
	    posxMod[i] = (int) PosXMod[i];
	    posyMod[i] = (int) PosYMod[i];  
      }
    }
  
  // Flags
  int FLAG_EXP_POFFTIME = 0; // indicates if the time at which the pacemaker stopped was exported.
  
  // Initial values
  double init_value_v	  = -1.0367; //Why these values? 
  double init_value_w 	  = -0.6656; //Why these values?
  double init_value_diff  = D;
  double init_value_ipace = 0.;
  double init_value_wl    = wL;
  double init_value_vt    = 0.;
  double init_value_wt    = 0.;
  
  // Arrays and their parameters
  mwSize num_rows = Ny;
  mwSize num_col  = Nx;
  
  typedef vector< vector<double> > TwoDVec;
  TwoDVec v;
  TwoDVec vt;
  TwoDVec w;
  TwoDVec wt; 
  TwoDVec diff;		// Diffusion
  TwoDVec Ipace;	// Pacemaker
  TwoDVec wl;		// w_L
  
  // Initialize the 2-D arrays
  v    .resize( num_col, vector<double>(num_rows, init_value_v)     );
  w    .resize( num_col, vector<double>(num_rows, init_value_w)     );
  diff .resize( num_col, vector<double>(num_rows, init_value_diff)  );
  /* Uncomment to set random diffusion coefficient at each pixel
  int iColumns, iRows;
  for(iColumns=0;iColumns<num_col;iColumns++){
    for(iRows=0;iRows<num_rows;iRows++){
	  diff[iColumns][iRows] = diff[iColumns][iRows] * (0.8 + (0.4*rand()/RAND_MAX));
	}
  }
  */
  Ipace.resize( num_col, vector<double>(num_rows, init_value_ipace) ); 
  wl   .resize( num_col, vector<double>(num_rows, init_value_wl)    );
  vt   .resize( num_col, vector<double>(num_rows, init_value_vt)    );
  wt   .resize( num_col, vector<double>(num_rows, init_value_wt)    );
  
  // Check convergence (CFL number for explicit scheme)
  //double nu = 2*D*dt/dx/dx;
  //if(nu > 0.5){
    //cout << "ERROR: divergent, adjust CFL number." << endl;
    //exit(0);
  //}  
  
//--------------------------------------------------------------//
    

  for(mwSize n=0; n < Nt; n++){	// LOOP OVER TIME
	  if(n<MapModTime){
		// Integration over space ------------------ // 
		#pragma omp parallel for
		{	
		  for(int k=0; k < nPosXY; ++k){
		  // Pacemakers
		  for(mwSize p=0; p<P.size(); p++){
			double dist_from_pace = sqrt( (P[p].x - posx[k])*(P[p].x - posx[k]) + (P[p].y - posy[k])*(P[p].y - posy[k]) ); //distance from pacemaker
			if(dist_from_pace <= P[p].r && n > P[p].on) //If the point is inside the pacemaker zone AND the pacemaker is turned on
			{	    
			  Ipace[posx[k]][posy[k]] = P[p].Ip; //Set the initial current in that pixel
			  wl[posx[k]][posy[k]] = P[p].wLP; //Set the inital wl in that pixel
			
			  if(n > P[p].off || CountBeats >= TotalBeats){ //If the pacemaker is off or the number of beats is reached
				Ipace[posx[k]][posy[k]] = 0; //Set current to zero
				wl[posx[k]][posy[k]] = wL; //Set wl to zero
			  }
			  
			  if(CountBeats == TotalBeats && FLAG_EXP_POFFTIME == 0){ //record time at which pacemaker is turned off
				outPace[0]  = n;
				FLAG_EXP_POFFTIME = 1;
			  }
			
			}
		  }	  
			
		// Model equations
		double AvgDiff = v[posx[k]+1][posy[k]] + v[posx[k]-1][posy[k]] + v[posx[k]][posy[k]+1] + v[posx[k]][posy[k]-1] - 4*v[posx[k]][posy[k]];
		double sigma   = (wh - wl[posx[k]][posy[k]])/(1 + exp(-4*v[posx[k]][posy[k]])) + wl[posx[k]][posy[k]];
		
		if(NoiseSignal==0){
			vt[posx[k]][posy[k]] = v[posx[k]][posy[k]] + dt * ( (v[posx[k]][posy[k]] - pow(v[posx[k]][posy[k]],3)/3 - w[posx[k]][posy[k]])/eps + Ipace[posx[k]][posy[k]]) + diff[posx[k]][posy[k]]*AvgDiff*dt/(dx*dx); 
			wt[posx[k]][posy[k]] = w[posx[k]][posy[k]] + dt * eps * sigma * (v[posx[k]][posy[k]] + beta - gamma*w[posx[k]][posy[k]]);
		}
		
		if(NoiseSignal==1){
			vt[posx[k]][posy[k]] = (0.98 + (0.04*rand()/RAND_MAX))*(v[posx[k]][posy[k]] + dt * ( (v[posx[k]][posy[k]] - pow(v[posx[k]][posy[k]],3)/3 - w[posx[k]][posy[k]])/eps + Ipace[posx[k]][posy[k]]) + diff[posx[k]][posy[k]]*AvgDiff*dt/(dx*dx)); 
			wt[posx[k]][posy[k]] = (0.98 + (0.04*rand()/RAND_MAX))*(w[posx[k]][posy[k]] + dt * eps * sigma * (v[posx[k]][posy[k]] + beta - gamma*w[posx[k]][posy[k]]));
		}
		
		// count beats at pacemaker
		if(posx[k]==P[0].x && posy[k]==P[0].y)
		{
		  if(vt[posx[k]][posy[k]] < threshold && v[posx[k]][posy[k]] >= threshold) // downwards
		  // If the new voltage is below threshold and the old is above
			CountBeats++; 
		}
		
		if(NoiseRecord==1){
		  if(n % SKIP_VDATA==0){
			outVid[posy[k] + posx[k]*Ny + Nx*Ny*n/SKIP_VDATA] = (0.6 + (0.8*rand()/RAND_MAX))*v[posx[k]][posy[k]]; 
		  }
		}
		
		else{
		  if(n % SKIP_VDATA==0){
			outVid[posy[k] + posx[k]*Ny + Nx*Ny*n/SKIP_VDATA] = v[posx[k]][posy[k]]; 
		  }
		}

		  } // loop over k
		} // loop pragma
		
		 #pragma omp parallel for
		{
		  // Update arrays
		  for(mwSize k=0; k < nPosXY; k++){
			v[posx[k]][posy[k]] = vt[posx[k]][posy[k]];
			w[posx[k]][posy[k]] = wt[posx[k]][posy[k]];
		  }
		}
	} //if(n < MapModTime)
	else{
		// Integration over space ------------------ // 
		#pragma omp parallel for
		{	
		  for(int k=0; k < nPosXYMod; ++k){
		  // Pacemakers
		  for(mwSize p=0; p<P.size(); p++){
			double dist_from_pace = sqrt( (P[p].x - posxMod[k])*(P[p].x - posxMod[k]) + (P[p].y - posyMod[k])*(P[p].y - posyMod[k]) ); //distance from pacemaker
			if(dist_from_pace <= P[p].r && n > P[p].on) //If the point is inside the pacemaker zone AND the pacemaker is turned on
			{	    
			  Ipace[posxMod[k]][posyMod[k]] = P[p].Ip; //Set the initial current in that pixel
			  wl[posxMod[k]][posyMod[k]] = P[p].wLP; //Set the inital wl in that pixel
			
			  if(n > P[p].off || CountBeats >= TotalBeats){ //If the pacemaker is off or the number of beats is reached
				Ipace[posxMod[k]][posyMod[k]] = 0; //Set current to zero
				wl[posxMod[k]][posyMod[k]] = wL; //Set wl to zero
			  }
			  
			  if(CountBeats == TotalBeats && FLAG_EXP_POFFTIME == 0){ //record time at which pacemaker is turned off
				outPace[0]  = n;
				FLAG_EXP_POFFTIME = 1;
			  }
			
			}
		  }	  
			
		// Model equations
		double AvgDiff = v[posxMod[k]+1][posyMod[k]] + v[posxMod[k]-1][posyMod[k]] + v[posxMod[k]][posyMod[k]+1] + v[posxMod[k]][posyMod[k]-1] - 4*v[posxMod[k]][posyMod[k]];
		double sigma   = (wh - wl[posxMod[k]][posyMod[k]])/(1 + exp(-4*v[posxMod[k]][posyMod[k]])) + wl[posxMod[k]][posyMod[k]];
		
		if(NoiseSignal==0){
			vt[posxMod[k]][posyMod[k]] = v[posxMod[k]][posyMod[k]] + dt * ( (v[posxMod[k]][posyMod[k]] - pow(v[posxMod[k]][posyMod[k]],3)/3 - w[posxMod[k]][posyMod[k]])/eps + Ipace[posxMod[k]][posyMod[k]]) + diff[posxMod[k]][posyMod[k]]*AvgDiff*dt/(dx*dx); 
			wt[posxMod[k]][posyMod[k]] = w[posxMod[k]][posyMod[k]] + dt * eps * sigma * (v[posxMod[k]][posyMod[k]] + beta - gamma*w[posxMod[k]][posyMod[k]]);
		}
		
		if(NoiseSignal==1){
			vt[posxMod[k]][posyMod[k]] = (0.98 + (0.04*rand()/RAND_MAX))*(v[posxMod[k]][posyMod[k]] + dt * ( (v[posxMod[k]][posyMod[k]] - pow(v[posxMod[k]][posyMod[k]],3)/3 - w[posxMod[k]][posyMod[k]])/eps + Ipace[posxMod[k]][posyMod[k]]) + diff[posxMod[k]][posyMod[k]]*AvgDiff*dt/(dx*dx)); 
			wt[posxMod[k]][posyMod[k]] = (0.98 + (0.04*rand()/RAND_MAX))*(w[posxMod[k]][posyMod[k]] + dt * eps * sigma * (v[posxMod[k]][posyMod[k]] + beta - gamma*w[posxMod[k]][posyMod[k]]));
		}
		
		// count beats at pacemaker
		if(posxMod[k]==P[0].x && posyMod[k]==P[0].y)
		{
		  if(vt[posxMod[k]][posyMod[k]] < threshold && v[posxMod[k]][posyMod[k]] >= threshold) // downwards
		  // If the new voltage is below threshold and the old is above
			CountBeats++; 
		}
		
		if(NoiseRecord==1){
		  if(n % SKIP_VDATA==0){
			outVid[posyMod[k] + posxMod[k]*Ny + Nx*Ny*n/SKIP_VDATA] = (0.6 + (0.8*rand()/RAND_MAX))*v[posxMod[k]][posyMod[k]]; 
		  }
		}
		
		else{
		  if(n % SKIP_VDATA==0){
			outVid[posyMod[k] + posxMod[k]*Ny + Nx*Ny*n/SKIP_VDATA] = v[posxMod[k]][posyMod[k]]; 
		  }
		}

		  } // loop over k
		} // loop pragma
		
		 #pragma omp parallel for
		{
		  // Update arrays
		  for(mwSize k=0; k < nPosXYMod; k++){
			v[posxMod[k]][posyMod[k]] = vt[posxMod[k]][posyMod[k]];
			w[posxMod[k]][posyMod[k]] = wt[posxMod[k]][posyMod[k]];
		  }
		}
		
	}//else
  
  } // loop over time
  
  //return 0;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
	/* initialize random seed: */
	srand (time(NULL)); //to call a random value between 0 and 1: rand()/RAND_MAX
    /* Inputs */
    mwSize Nt; 
    mwSize pon, poff;
    mwSize pdx, pdy; //position to pacemaker, relative to center
    mwSize TotalBeats;
    mwSize nPosXY, nPosXYMod;
    double D, dx, dt;
    double wH, wLP;
    double *PosX, *PosY, *PosXMod, *PosYMod;
    int MapModTime;
    int NoiseSignal;
    int NoiseRecord;
    
    /* output matrix */
    double *outVid, *outPace;         

    /* check for proper number of arguments */
    if(nrhs!=20) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","20 inputs required.");
    }

    /* make sure the first input argument is scalar */
    for(int i=0; i<11; i++)
    {
      if( !mxIsDouble(prhs[i]) || 
         mxIsComplex(prhs[i]) ||
         mxGetNumberOfElements(prhs[i])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Inputs must be a scalars.");
      }
    }
    
    
        
    /* get the inputs  */
    Nt 		= mxGetScalar(prhs[0]);
    dx 		= mxGetScalar(prhs[1]);
    dt 		= mxGetScalar(prhs[2]);
    D 		= mxGetScalar(prhs[3]);
    wH	 	= mxGetScalar(prhs[4]);
    wLP 	= mxGetScalar(prhs[5]);
    pon 	= mxGetScalar(prhs[6]);
    poff 	= mxGetScalar(prhs[7]);
    pdx 	= mxGetScalar(prhs[8]);
    pdy 	= mxGetScalar(prhs[9]);
    TotalBeats 	= mxGetScalar(prhs[10]);
    nPosXY 	= mxGetScalar(prhs[11]);
    PosX	= mxGetPr(prhs[12]); 
    PosY	= mxGetPr(prhs[13]); 
    NoiseSignal = mxGetScalar(prhs[14]);
    NoiseRecord = mxGetScalar(prhs[15]);
    nPosXYMod 	= mxGetScalar(prhs[16]);
    PosXMod	= mxGetPr(prhs[17]); 
    PosYMod	= mxGetPr(prhs[18]); 
    MapModTime = mxGetScalar(prhs[19]);
    
    /* set the output pointer to the output matrix */
    mwSize dimVid[3] = {82, 82, Nt/SKIP_VDATA};
    
    plhs[0] = mxCreateNumericArray(3, dimVid, mxDOUBLE_CLASS, mxREAL);	// video
    plhs[1] = mxCreateDoubleMatrix(10, 1, mxREAL); 			// info
    
    /* create a C pointer to a copy of the output matrix */
    outVid    = mxGetPr(plhs[0]);
    outPace   = mxGetPr(plhs[1]);
    
    /* call the computational routine */
    myspiral(Nt, dx, dt, D, wH, wLP, pon, poff, pdx, pdy, TotalBeats, 
		    nPosXY,
		    PosX, PosY, NoiseSignal, NoiseRecord, nPosXYMod, PosXMod, PosYMod, MapModTime,
		    outVid, outPace); 
}
