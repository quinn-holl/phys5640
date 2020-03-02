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
#include "TRandom3.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TString.h"

const int N = 4; //Number of spheres in simulation
const int TIME = 5000; //How long system will run in # of collisions
const int SUBTIME = 5; //number of updates between collisions
//const double sigma = 0.075; //Radius of hard spheres
const double INF = 1000000; //Definition of Infinity here to compare times
const int DIM = 2; //we're only doing 2d Here
const double SIZE = 1; //Box is of length size
const double VABS = .1; //This is just part of the initial condition
const int THERMAL = 1000; //after how many events do we start sampling
const double LIMIT = 0.15; //Largest size of sigma
const int NEXP = 5;
const int NTRIALS = 5000;

struct System{
  std::array<double,N> x;
  std::array<double,N> y;
  std::array<double,N> vx;
  std::array<double,N> vy;
}; 

std::array<double,DIM> generatePoint(std::array<double,N> x,std::array<double,N> y,double sigma,TRandom *r3,int BCflag,int index);
std::array<double,DIM> box(std::array<double,DIM> x);
std::array<double,DIM> diffV(std::array<double,DIM> x,std::array<double,DIM> y);
std::vector<double> plots(double sigma,int BCflag,TRandom *r3);


int main(int argc, char **argv){
  TApplication theApp("App",&argc,argv);
  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  /**  Outline for Ex 2.6
       Run NEXP, or the number of different covering densities. For each density
       Run NTRIALS and for each generate 4 legal points of the centers for the
       disks. Add each generated disk's coordinate to the histogram. 
  */


  

  
  std::array<double,NEXP> sig = {0.01,0.025,0.05,0.1,0.15}; //sigma values
  TRandom *r3 = new TRandom();
  std::array<std::vector<double>,NEXP> hists;
  std::array<std::vector<double>,NEXP> histsPeriod;
  for(int j = 0;j<NEXP; j++){
    hists[j] = plots(sig[j],0,r3); //box boundary conditions
    histsPeriod[j] = plots(sig[j],1,r3); //Periodic BC
  }

  TCanvas *cr[NEXP];
  for(int i = 0;i<NEXP;i++){
    TString s;
    s.Form("c%d",i);
    cr[i] = new TCanvas(s,"Canvas",dw,dh);
    cr[i]->Divide(2,1);
    cr[i]->cd(1);
    TString s1,s2;
    s1.Form("Box BC Direct Sample #sigma = %lf;X-val;Count",sig[i]);
    s2.Form("Periodic BC Direct Sample #sigma = %lf;X-val;Count",sig[i]);
    TH1F *h1 = new TH1F("h1",s1,100,0.0,1.0);
    TH1F *h2 = new TH1F("h2",s2,100,0.0,1.0);
    printf("SIGMA: %e \n",sig[i]);
    for(int j = 0;j<NTRIALS*N;j++){
      printf("%d\n",j);
      double value1 = hists[i][j];
      double value2 = histsPeriod[i][j];
      h1->Fill(value1);
      h2->Fill(value2);
    }
    h1->Draw();
    cr[i]->cd(2);
    h2->Draw();
  }
  
  
  



  std::cout << "Press ^c to exit" << std::endl;
  theApp.SetIdleTimer(120,".q");
  theApp.Run();
}

std::array<double,DIM> generatePoint(std::array<double,N> x,std::array<double,N> y,double sigma,TRandom *r3,int BCflag,int index){
  double xtemp = 0.0;  //If debugging, check that were not just getting 
  double ytemp = 0.0; // disks at 0,0
 
  if(BCflag == 0){
    xtemp = r3->Rndm()*(SIZE-2*sigma) + sigma;
    ytemp = r3->Rndm()*(SIZE-2*sigma) + sigma;
  } else {
    xtemp = r3->Rndm()*SIZE;
    ytemp = r3->Rndm()*SIZE;
  }
  if(index == -1){
    std::array<double,DIM> p1 = {xtemp,ytemp};
    return p1;
  }
  for(int i = 0;i<=index;i++){ //could make speed improvements here
    double distance = 0.0;
    if(BCflag == 0){
      distance = pow(pow(xtemp-x[i],2) + pow(ytemp-y[i],2),0.5);
    }
    if(BCflag == 1){
      std::array<double,DIM> diffVector;
      std::array<double,DIM> p1 = {x[i],y[i]};
      std::array<double,DIM> p2 = {xtemp,ytemp};
      diffVector = diffV(p1,p2);
      distance = pow(pow(diffVector[0],2) + pow(diffVector[1],2),0.5);
    }
    if( fabs(distance -1E-5)< 1E-5){
      printf("Something went wrong using distance zero\n");
    }
    if(distance < 2*sigma){
      std::array<double,DIM> retval = {-1,-1};
      return retval;
    }
  }
  std::array<double,DIM> retval = {xtemp,ytemp};
  return retval;
}

std::array<double,DIM> box(std::array<double,DIM> x){
  std::array<double,DIM> retval;
  retval[0] = std::fmod(x[0],SIZE);
  retval[1] = std::fmod(x[1],SIZE);
  if(retval[0] < 0){
    retval[0] = retval[0] + SIZE;
  } 
  if(retval[1] < 1){
    retval[1] = retval[1] + SIZE;
  }
  return retval;
}

std::array<double,DIM> diffV(std::array<double,DIM> x,std::array<double,DIM> y){
  std::array<double,DIM> del;
  del[0] = y[0]-x[0];
  del[1] = y[1]-x[1];
  del = box(del);
  if(del[0] > SIZE/2.0){
    del[0] = del[0] - SIZE/2.0;
  }
  if(del[1] > SIZE/2.0){
    del[1] = del[1] -SIZE/2.0;
  }
  return del;
}


std::vector<double> plots(double sigma,int BCflag,TRandom *r3){
  std::vector<double> sigHist;
  int l = 0;
  while(l<NTRIALS){
    int rejFlag = 0;
    std::array<double,N> x;
    std::array<double,N> y;
    for(int i = 0; i<N;i++){ //generate 4 ind points
      std::array<double,DIM> p1 = generatePoint(x,y,sigma,r3,BCflag,i-1);
      if( fabs(p1[0] - (-1)) < 1E-5){
	rejFlag = 1;
	break;
      }
      x[i] = p1[0];
      y[i] = p1[1];
    }
    if(rejFlag == 0){ //essentially if our config was accepted
      for(int m = 0;m<N;m++){
	sigHist.push_back(x[m]);
      }
      l++; //now we can increment another trial;
    } 
  }
  return sigHist;
}
