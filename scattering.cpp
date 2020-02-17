#include <cmath>
#include "TGraph.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <array>
#include "TLegend.h"
#include "TMultiGraph.h"

const double EPSILON = 5.9; //milli eV, meV
const double RHO = 3.57;  //Angstrom
const double K = 6.12; //(meV)^-1 (rho)^-2  and is equivalent to 2m/hbar^2
const double C = 1; //HO example where h^2/m = 1 
const int N = 1000; //How many points in the graph for LJ potential
#define STEPS 1000


using namespace std;

double ljPotential(double radius);
double firstStep(double psi0,double psi0prime,double h,double Energy,double r,
		 double l);
std::array<std::array<double,2>,STEPS> Numerov(double psi0,double psi0prime,
					       double steps, double a,double b,
					       double Energy, double l);
double kVector(double r,double Energy,double l);


int main(int argc, char **argv){
  //Want to plot the lj potential in main over the range 3-10 (Angstrom)
  TApplication theApp("App",&argc,argv);
  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  TCanvas *c1 = new TCanvas("c1","Canvas",dw,dh);
  TGraph *g = new TGraph();
  TGraph *expected = new TGraph();

  double psi0 = 0;
  double psi1= pow(r,l+1);
  double a = 0.0;
  double b = 4.0;
  std::array<std::array<double,2>,STEPS> sol = Numerov(psi0,psi1,STEPS,a,b,3.0,0);
  double h = (a-b)/STEPS;
  double Interror = 0.0;
  for(unsigned int i = 0;i<sol.size();i++){
    g->SetPoint(g->GetN(),sol[i][0],sol[i][1]);
    double r = sol[i][0];
    double expectedVal = r*exp(h*h/2)*exp(-r*r/2);
    Interror += fabs(expectedVal - sol[i][1]);
    expected->SetPoint(expected->GetN(),r,expectedVal);
  }
  printf("Total Integrated Error For HO WF: %e\n",Interror);
  c1->cd(1);
  g->SetLineColor(2);
  expected->SetLineColor(3);
  g->SetTitle("Radial Wavefunction 3d-HO model;radius;psi");
  g->Draw("a*");
  expected->SetTitle("Analytic Solution");
  expected->Draw("same");
  c1->BuildLegend();

  /*  If you want to create the LJ plot, uncommment out
  TGraph *g = new TGraph();
  double start = 3;  //Units are in Angstrom
  double end = 10;   //Units are in Angstrom
  for(int i = 0;i<N;i++){
    double x = start +(end-start)*((double)i)/N;
    double y = ljPotential(x);
    g->SetPoint(g->GetN(),x,y);
  }

  g->SetTitle("LJ Potential; Radius (Angstrom); (meV)");
  g->SetLineColor(2);
  c1->cd(1);
  g->Draw("ca");
  */

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(60,".q");
  theApp.Run();
  
}


double ljPotential(double radius){
  /**  Uncomment out for Lennard Jones Potential
   */
  double term1 = pow((RHO/radius),6);
  double term2 = pow(term1,2);
  double retval = EPSILON*(term2-2*term1);
  //double retval = pow(radius,2);
  return retval;
}

double firstStep(double psi0,double psi0prime,double h,double Energy,double r,
		 double l){
  double term1 = psi0*(1-(1.0/24.0)*(h*h*kVector(r+2*h,Energy,l)));
  double term2 = h*psi0prime*(1-(1.0/12.0)*(h*h*kVector(r+2*h,Energy,l)));
  double term3 = (h*h)/(24.0)*(7*kVector(r,Energy,l)*psi0);
  double term4 = -(pow(h,4))*kVector(r+2*h,Energy,l)/36.0*(kVector(r,Energy,l)*psi0);
  double term5 = (1-((h*h/4.0)*kVector(r+h,Energy,l))+((pow(h,4)/18.0)*kVector(r+h,Energy,l)*kVector(r+2*h,Energy,l)));
  double psi1 = (term1 + term2 + term3 -term4)/ term5;
  return psi1;
}

std::array<std::array<double,2>,STEPS> Numerov(double psi0,double psi1,
					       double steps, double a,double b,
					       double Energy,double l){

  //need to add l, in order to get magnitude of centrifugal force
  std::array<std::array<double,2>,STEPS> values;
  double numerov = 0.0;
  double h = (b-a)/steps;

  //First Step Calculation
  double r = a;
  

  for(int i = 0; i < steps; i++){
    r = (a+2*h) + i*h;
    double k1 = kVector(r-h,Energy,l);
    double k2 = kVector(r-2*h,Energy,l);
    double k3 = kVector(r,Energy,l);
    double numerator = 2*(1-(5.0/12.0)*h*h*k1)*psi1 -(1+(1.0/12.0)*h*h*k2)*psi0;
    double denominator = 1 + (1.0/12.0)*h*h*k3;    
    numerov = (numerator/denominator);
    values[i][0] = r;
    values[i][1] = numerov;
    psi0 = psi1;
    psi1 = numerov;
  }

  return values;
}


double kVector(double r,double Energy,double l){
  if(r==0){
    return 0;
  }
  double integrand = C*(Energy - ljPotential(r) - (l*(l+1))/(r*r));
  printf("%e    %e  \n",Energy,ljPotential(r));
  printf("%e\n",integrand);
  return integrand;
}
