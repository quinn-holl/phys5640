// Calculate the Voltage for a 2D system with wire at fixed voltage immediately above
// a 3-sided box at 0V

#include "TGraph2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBox.h"
#include "TApplication.h"

#include <getopt.h>
#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;

const double XBEGIN = 0.35;
const double XEND = 0.65;
const double YBEGIN = 0.45;
const double YEND = 0.55;
const double RHO = 3.33;
const double RSTART = 0.15;
const double REND = 0.25;
const double PHISTART = 0; //degrees to radians
const double PHIEND = 35*0.0175; 
const int THICK = 1; //leads to 2*THICK thickness of capacitor plates



// generic code to do one iteration of finite difference method
// Jacobi Method



double iterateJ(vector<vector<double>> &V, double L,int BC){
  int xstart = int(V.size()*XBEGIN);
  int xend = int(V.size()*XEND);
  int ystart = int(V.size()*YBEGIN);
  int yend = int(V.size()*YEND);
  auto Vtmp = V;
  double dVmax=1e-50;
  int nx=V.size();
  int ny=V[0].size();
  double delta = L/(V.size()-1);

  int thick = 0;
  if(BC == 3){
    thick =1 ;
  }

  for (int i=1; i<nx-1; i++){
    for (int j=1; j<ny-1; j++){
      double Vnew = 0.25*(Vtmp[i+1][j]+Vtmp[i-1][j]+Vtmp[i][j+1]+Vtmp[i][j-1]);
      double dV=fabs(Vnew-V[i][j]);
      if( (i >= xstart && i <= xend) && (  (j>= ystart-thick) && (j <= ystart+thick))){
	if(BC == 2){
	  Vnew += M_PI*RHO*delta;
	  V[i][j] = Vnew;
	}
	continue;
      }
      if((i >= xstart && i <= xend) && ( (j >= yend-thick) && (j<= yend+thick))){
	if(BC ==2){
	  Vnew -= M_PI*RHO*delta;
	  V[i][j] = Vnew;
	}
	continue;
      }
      else{
	dVmax=std::max(dVmax,dV);    // keep track of max change in this sweep
	V[i][j] = Vnew;
      }
    }
  }
  return dVmax;
}




//Uncomment out if you want cylindircal solver 


double iterateJCyl(vector<vector<double>> &V, double L){
  if(V.size() == 0 || V[0].size() == 0){
    printf("Size of array is off\n");
    exit(1);
  }
  double deltaR = 2.0/(2.0*V.size()+1);
  double deltaP = 2.0*M_PI/(V[0].size());
  int rstart = int(RSTART/deltaR);
  int rend = int(REND/deltaR);
  int phiStart = int(PHISTART/deltaP);
  int phiEnd = int(PHIEND/deltaP);
  auto Vtmp = V;
  double dVmax = 1e-50;
  int nr = V.size();
  int np = V[0].size();
  for(int i = 1; i<nr-1;i++){
    for(int j = 0;j < np; j++){
      double r = (i-0.5)*deltaR;
      double V1 =  (V[i+1][j] + V[i-1][j])/(deltaR*deltaR);
      double V2 = (V[i+1][j] - V[i-1][j])/(2*deltaR*r);
      // Heed the BC for cylindrical, that theta wraps back around itself
      double V3 = (V[i][j+1] + V[i][j-1])/(deltaP*deltaP*r*r);
      if(j == 0){
	V3 = (V[i][j+1] + V[i][np-1])/(deltaP*deltaP*r*r);
      } else if (j == np-1){
	V3 = (V[i][0] + V[i][j])/(deltaP*deltaP*r*r);
      }
      double Vnew = (V1+V2+V3)/((2/(pow(deltaR,2))) + (2/(pow(deltaP*r,2))));
      double dV = fabs(Vnew-V[i][j]);
      if( (j>= phiEnd || j <= phiStart) && (i== rstart || i ==rend)){
	continue; //keeping constant potential
      }
      else if( j == np-1 ){  //B.C. for cylindrical coordinates
	V[i][j]  = V[i][0];
      }
      else{
	dVmax = std::max(dVmax,dV);
	V[i][j] = Vnew;
      }
    }
  }
  return dVmax;


}

void printLattice(vector<vector<double>> &V){
  for(int i = 0; i<V.size();i++){
    for(int j = 0;j<V[0].size();j++){
      printf("%.2lf \t ",V[i][j]);
    }
    printf("\n");
  }
}



// Gauss-Seidel Method
double iterateGS(vector<vector<double>> &V){
  double dVmax=1e-50;
  int nx=V.size();
  int ny=V[0].size();
  for (int i=1; i<nx-1; i++){
    for (int j=1; j<ny-1; j++){
      double Vnew = 0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]);
      double dV=fabs(Vnew-V[i][j]);
      dVmax=std::max(dVmax,dV);    // keep track of max change in this sweep
      V[i][j] = Vnew;
    }
  }
  return dVmax;
}

// fill a TGraph2D object from a vector of voltages
// delta: grid spacing
// the optional range parameter defines the subregion to plot

void fillGraph(TGraph2D* tg, const vector<vector<double>> &V, double delta, TBox *range=0){
  int nx=V.size();
  int ny=V[0].size();
  tg->Clear();                 // reset the graph
  for (int i=0; i<nx; i++){
    double x = i*delta;
    for (int j=1; j<ny; j++){
      double y = j*delta;
      tg->SetPoint(tg->GetN(),x,y,V[i][j]);
    }
  }
}

// Define box 0<x<L, 0<y<L
// eps: convergence criteria (max size of change at any grid point in an iteration)
// maxIter: max iterations in case of non-converence
// Npts : smoothness parameter, number of grid points in x,y
// pass a tcanvas for an animated solution, with specified max rate of frames/second
TGraph2D* LaplaceLine(int maxIter=500, double eps=0.001, int Npts=100, int Mpts = 100, TCanvas *tc=0, int rate=10, int BC = 1){
  double L=100;            // length of any side
  double Vtop= 100;         // Voltage at top of box
  double Vbot = -100;
  int maxgraphlines=200;   // max lines to draw in each direction
  vector<vector<double>> V(Npts, vector<double> (Npts, 0));  // create N x N vector, init to 0
  int xstart = int(V.size()*XBEGIN);
  int xend = int(V.size()*XEND);
  int ystart = int(V.size()*YBEGIN);
  int yend = int(V.size()*YEND);
  double delta = L/(Npts-1);                                 // grid spacing

  int thick = 0;
  if(BC == 3){
    thick = 1;
  }
  for(int i = xstart; i<= xend; i++){
    for(int j = -1*thick; j<=thick; j++){
      if(BC == 1 || BC == 2 || BC == 3){
	V[i][ystart+j] = Vtop;
	V[i][yend+j] = Vbot;
      }
      if(BC == 4){
	double w = xend - xstart;
	double x = i - xstart;
	if(x <= 0.5*w){
	  V[i][ystart+j] = 200*x/w;
	} else {
	  V[i][ystart+j] = 200*(1-x/w);
	}
	V[i][yend + j] = -V[i][ystart +j];
      }
      if(BC == 5){
	double w = xend-xstart;
	double x = i -xstart;
	V[i][ystart+j] = 100*sin(2*M_PI*x/w);
	V[i][yend + j]= -V[i][ystart+j];
      }
    }
  }
  int msec = 1000/rate;                                      // milliseconds sleep between frames
  TBox *plotRange = new TBox(0,0,1.1*L,1.1*L);

  TGraph2D* tgV = new TGraph2D();                            // graph to store result
  if (Npts<50) tgV->SetLineWidth(3);                         
  tgV->SetLineColor(kRed);
  tgV->SetNpx(std::min(maxgraphlines,Npts));  tgV->SetNpy(std::min(maxgraphlines,Npts)); 
  tgV->SetTitle("Voltage;x;y;V");
  
  double dV;
  int niter=0;
  do{
    dV=iterateJ(V,L,BC);   // iterate using Jacobi method
    //dV=iterateGS(V);   // iterate using Gauss-Seidel method
    ++niter;
    /**
    if (tc) {
      tc->cd();
      fillGraph(tgV,V,delta,plotRange);
      tgV->Draw("surf");
      tc->Update();
      gSystem->Sleep(msec);
    }
    */
  } while (dV>eps && niter<maxIter);

  //Q3, create a graph of the charge distributions of the capacitors
  if(BC == 3){
    tgV->Clear();
    for(int i = 0; i<(xend-xstart+1);i++){
      for(int j = 0;j<(2*THICK+1);j++){
	int indexX = xstart+i;
	int indexY = ystart-THICK+j;
	double Ptop = (V[indexX][indexY]-0.25*(V[indexX+1][indexY] + V[indexX-1][indexY] + V[indexX][indexY+1] + V[indexX][indexY-1]))/(M_PI*delta*delta);
	tgV->SetPoint(tgV->GetN(),indexX,indexY,Ptop);
      }
    }
  }
  
  cout << "Ended calculation with " << niter << " iterations, dVmax = " << dV << endl;
  if(BC != 3){
    fillGraph(tgV,V,delta,plotRange);  
  }
  return tgV;
}



//Uncomment out if you want cylindrical solver

void fillGraphCyl(TGraph2D* tg, const vector<vector<double>> &V, double deltaR, double deltaP, TBox *range=0){
  int nx=V.size();
  int ny=V[0].size();
  tg->Clear();                 // reset the graph
  for (int i=1; i<nx; i++){
    double r = (i-0.5)*deltaR;
    for (int j=1; j<ny; j++){
      double theta = (j-1)*deltaP;
      double x = r*cos(theta);
      double y = r*sin(theta);
      tg->SetPoint(tg->GetN(),x,y,V[i][j]);
    }
  }
}

TGraph2D* LaplaceLineCyl(int maxIter = 100, double eps = 0.001, int Npts = 100, int Mpts = 100, TCanvas *tc =0, int rate = 10, int BC = 6){
  double L = 0;
  double Vin = -100;
  double Vout = 100;
  int maxgraphlines = 200;
  vector<vector<double>> V(Npts,vector<double> (Mpts,0));
  double deltaR = 2.0/(2.0*V.size()+1);
  double deltaP = 2.0*M_PI/V[0].size();
  int rstart = int(RSTART/deltaR);
  int rend = int(REND/deltaR);
  int phiStart = int(PHISTART/deltaP);
  int phiEnd = int(PHIEND/deltaP);
  for( int i = 0; i<Mpts;i++){
    V[Npts-1][i] = 0;//creating that box for grounded B.C. conditions
  }
  for(int j = 0; j<Mpts; j++){
    if( j>= phiStart && j <= phiEnd){
      continue;
    } else {
      V[rstart][j] = Vin;
      V[rend][j] = Vout;
    }
  }
  int msec = 1000/rate;
  TBox *plotRange = new TBox(-1.1*L,-1.1*L,1.1*L,1.1*L);
  TGraph2D *tgV = new TGraph2D();
  if (Npts<50) tgV->SetLineWidth(3);                         
  tgV->SetLineColor(kRed);
  tgV->SetNpx(std::min(maxgraphlines,Npts)); 
  tgV->SetNpy(std::min(maxgraphlines,Npts)); 
  tgV->SetTitle("Voltage;x;y;V");
  
  double dV;
  int niter=0;
  do{
    dV=iterateJCyl(V,L);   // iterate using Jacobi method
    ++niter;
    if (tc) {
      tc->cd();
      fillGraphCyl(tgV,V,deltaR, deltaP,plotRange);
      tgV->Draw("surf");
      tc->Update();
      gSystem->Sleep(msec);
    }
  } while (dV>eps && niter<maxIter);
  
  cout << "Ended calculation with " << niter << " iterations, dVmax = " << dV << endl;
  //Draw the Electric Field Strength Plot, again, only showing magnitude not direction
  TGraph2D *tgE = new TGraph2D();
  tgE->SetTitle("Electric Field Magnitude");
  TCanvas *c3 = new TCanvas("c3","c3",1);
  c3->cd(1);
  deltaR = 2.0/(2.0*V.size()+1);
  deltaP = 2.0*M_PI/(V[0].size());
  int nr = V.size();
  int np = V[0].size();
  vector<vector<double>> E(Npts,vector<double> (Mpts,0));
  for(int i = 1; i<nr-1;i++){
    for(int j = 0;j < np; j++){
      double r = (i-0.5)*deltaR;
      double E1 = (V[i+1][j]-V[i-1][j])/(2*deltaR);
      // Heed the BC for cylindrical, that theta wraps back around itself
      double E2 = (V[i][j+1]-V[i][j-1])/(2*r*deltaP);
      if(j == 0){
	E2 = (V[i][j+1] - V[i][np-1])/(2*r*deltaP);
      }else if (j == np-1){
	E2 = (V[i][0] - V[i][j-1])/(2*r*deltaP);
      }
      double Emag = pow(  pow(E1,2) + pow(E2,2) , 0.5);
      if( j == np-1 ){  //B.C. for cylindrical coordinates
	E[i][j]  = E[i][0];
      }
      else{
	E[i][j] = Emag;
      }
    }
  }
  fillGraphCyl(tgE,E,deltaR,deltaP,plotRange);
  tgE->Draw("SURF3");

 
  fillGraphCyl(tgV,V,deltaR, deltaP,plotRange);
  return tgV;
}









void usage(char *prog){
  std::cerr << "Usage: " << prog << " <option(s)> SOURCES"
	    << "Options:\n"
	    << "\t-h\t\tShow this help message\n"
	    << "\t-a\t\tDisplay animation of solution"
    	    << "\t-I\t\t(max) Number of iterations [100]"
	    << "\t-e\t\tconvergence criteria [0.001]"
    	    << "\t-N\t\tNumber of points in x,y [100]"
	    << "\t-R\t\tmax frames/second with animation [10]"
            << "\t-p\t\tproblem # associated with wanted output"
	    << std::endl;
  exit(0);
}


int main(int argc, char *argv[]){
  TApplication theApp("App", &argc, argv, NULL, -1);  // -1 disables ROOT arg processing

  // defaults for LaplaceLine
  int maxIter=100;
  double eps=0.001;
  int Npts=100;
  TCanvas *tc=0;
  int rate=10;
  int BC = 1;
  
  int opt;
  while ((opt = getopt(argc, argv, "haI:e:N:r:p:")) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      break;
    case 'a':
      tc=new TCanvas();
      break;
    case 'I':
      maxIter=atoi(optarg);
      break;
     case 'e':
      eps=atof(optarg);
      break; 
    case 'N':
      Npts=atoi(optarg);
      break;
    case 'r':
      rate=atoi(optarg);
      break;
    case 'p':
      BC = atoi(optarg);
    }
  }

  printf("BC seen is %d  \n",BC);

  if(BC != 6){
    auto tg =LaplaceLine(maxIter,eps,Npts,Npts,tc,rate,BC);
    if(!tc) tc = new TCanvas();
    if(BC == 3){
      tg->Draw("COLZ");
    } else {
      tg->Draw("SURF3");
    }
  } else {
    auto tg = LaplaceLineCyl(maxIter,eps,Npts,Npts,tc,rate,BC);
    if(!tc) tc = new TCanvas();
    tg->Draw("SURF3");
  }
  if(BC == 4){
    TCanvas *c2 = new TCanvas("c2","c2",1);
    c2->cd(1);
    auto tg2 = LaplaceLine(maxIter,eps,Npts,Npts,c2,rate,5); //do the second periodic BC
    tg2->Draw("SURF3");
  }
  
               // explore other drawing options!
  //SURF3 for every question except 3
  //COlZ for question 3

  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
