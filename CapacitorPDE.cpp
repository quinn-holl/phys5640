// Calculate the Voltage for a 2D system with two capacitor plates

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

// generic code to do one iteration of finite difference method
// Jacobi Method
double iterateJ(vector<vector<double>> &V){
  int xstart = int(V.size()*XBEGIN);
  int xend = int(V.size()*XEND);
  int ystart = int(V.size()*YBEGIN);
  int yend = int(V.size()*YEND);
  auto Vtmp = V;
  double dVmax=1e-50;
  int nx=V.size();
  int ny=V[0].size();
  for (int i=1; i<nx-1; i++){
    for (int j=1; j<ny-1; j++){
      double Vnew = 0.25*(Vtmp[i+1][j]+Vtmp[i-1][j]+Vtmp[i][j+1]+Vtmp[i][j-1]);
      double dV=fabs(Vnew-V[i][j]);
      dVmax=std::max(dVmax,dV);    // keep track of max change in this sweep
      if( (i >= xstart && i <= xend) && ( j == ystart || j == yend)){
	continue;   //dont alter voltage of capacitors
      }
      else{
	V[i][j] = Vnew;
      }
    }
  }
  return dVmax;
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
      if (range && range->IsInside(x,y))
	tg->SetPoint(tg->GetN(),x,y,V[i][j]);
    }
  }
}

// Define box 0<x<L, 0<y<L
// eps: convergence criteria (max size of change at any grid point in an iteration)
// maxIter: max iterations in case of non-converence
// Npts : smoothness parameter, number of grid points in x,y
// pass a tcanvas for an animated solution, with specified max rate of frames/second
TGraph2D* LaplaceLine(int maxIter=100, double eps=0.001, int Npts=100, TCanvas *tc=0, int rate=10){
  double L=100;            // length of any side
  double Vtop=100;         // Voltage at top of box
  double Vbot = -100;
  int maxgraphlines=200;   // max lines to draw in each direction
  vector<vector<double>> V(Npts, vector<double> (Npts, 0));  // create N x N vector, init to 0
  int xstart = int(V.size()*XBEGIN);
  int xend = int(V.size()*XEND);
  int ystart = int(V.size()*YBEGIN);
  int yend = int(V.size()*YEND);
  double delta = L/(Npts-1);                                 // grid spacing
  for (int i=0; i<Npts; i++) {   //B.C for parallel plate capacitor
    V[i][Npts-1] = 0;
    V[i][0] = 0;
    V[Npts-1][i] = 0;
    V[0][i] = 0;
  }
  for(int i = xstart; i<= xend; i++){
    V[i][ystart] = Vtop;
    V[i][yend] = Vbot;
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
    dV=iterateJ(V);   // iterate using Jacobi method
    //dV=iterateGS(V);   // iterate using Gauss-Seidel method
    ++niter;
    if (tc) {
      tc->cd();
      fillGraph(tgV,V,delta,plotRange);
      tgV->Draw("surf");
      tc->Update();
      gSystem->Sleep(msec);
    }
  } while (dV>eps && niter<maxIter);
  
  cout << "Ended calculation with " << niter << " iterations, dVmax = " << dV << endl;
 
  fillGraph(tgV,V,delta,plotRange);
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
  
  int opt;
  while ((opt = getopt(argc, argv, "haI:e:N:r:")) != -1) {
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
    }
  }
 
  auto tg=LaplaceLine(maxIter,eps,Npts,tc,rate);

  // display final result
  if (!tc) tc=new TCanvas();
  tg->Draw("surf");              // explore other drawing options!
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
