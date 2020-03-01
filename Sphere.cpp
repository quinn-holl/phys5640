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
#include "TProfile.h"

const int N = 10; //Number of spheres in simulation
const int TIME = 10; //How long system will run in # of collisions
const int SUBTIME = 5; //number of updates between collisions
const double sigma = 0.075; //Radius of hard spheres
const double INF = 1000000; //Definition of Infinity here to compare times
const int DIM = 2; //we're only doing 2d Here
const double SIZE = 1; //Box is of length size
const double VABS = .1; //This is just part of the initial condition
const int THERMAL = 1; //after how many events do we start sampling

//Box Dimensions: [0,1],[0,1] for both x and y


//This struct should contain the entire state of our system at time t
//It has four arrays of size N, so how many particles are in our simulation
//Each array contains either position or velocity information, in one of 2 
//coordinates. Each sphere's state can be accessed in the following way
//One sphere is entirely specified by x[i],y[i],vx[i],vy[i]. All you need to 
//make sure is that you input index consistently.
struct System{
  std::array<double,N> x;
  std::array<double,N> y;
  std::array<double,N> vx;
  std::array<double,N> vy;
}; 

//outputs time till next pair coll
std::vector<double> pairCT(System states);  
System pairC(System states,int index1,int index2);  //updates system after a pair
std::vector<double>  wallCT(System states);  //time till next wall c
System wallC(System states,int index,int direction);   //updates system after a wall
System update(System states,double t); //propogates system through time t
double wallTime(double x, double vx); //one-coordinate calculation
System Initialize();
double pairTime(std::array<double,DIM> pairx,std::array<double,DIM> pairy,
		std::array<double,DIM> pairvx,std::array<double,DIM> pairvy);
std::array<double,DIM> generatePoint(std::array<double,N> x,std::array<double,N> y);





int main(int argc, char **argv){
  TApplication theApp("App",&argc,argv);
  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  TCanvas *c1 = new TCanvas("c1","Canvas",dw,dh);
  double t = 0; //initialize time variable
  System states = Initialize(); //Need to initialize N hard disks
  TH2F *hist = new TH2F("hist","Event Disk Hist",100,0.0,1.0,100,0.0,1.0);
 

  for(int i = 0;i<TIME;i++){  //this is in "collision time"
    std::vector<double> pairVector = pairCT(states);  
    double tpair = pairVector[0];
    int tpairIndex1 = std::round(pairVector[1]); //need to know which two disks
    int tpairIndex2 = std::round(pairVector[2]); // to update
    std::vector<double> wallVector = wallCT(states);
    double twall = wallVector[0];  
    int twallIndex = std::round(wallVector[1]);  //same goes here, need to
    int direction = std::round(wallVector[2]);
    double tnext = std::min(tpair,twall);  //know which disk to update
    if(tnext == INF){
      printf("Something's Gone Wrong, Using INF for time step\n");
    }
    for(int m = 1;m<= SUBTIME; m++){
      double deltaT = (1.0/SUBTIME)*(tnext);
      states = update(states,deltaT);
      for(int l = 0;l<N;l++){
	printf("collision: %d, particle: %d, x: %e  y:%e  time: %e\n",i,l,states.x[l],states.y[l],deltaT);
      }
      //probably want to sample here for histogram points
      if(i>THERMAL){
	for(int j = 0;j<N;j++){
	  hist->Fill(states.x[j],states.y[j]);
	}
      }
    }
    if(twall > tpair){ //updates states to reflect the collision at tnext
      states = pairC(states,tpairIndex1,tpairIndex2);
    } else {
      printf("vx:%e  vy:%e  x:%e   y:%e\n",states.vx[twallIndex],states.vy[twallIndex],states.x[twallIndex],states.y[twallIndex]);
      states = wallC(states,twallIndex,direction);
      printf("vx:%e  vy:%e  x:%e   y:%e\n",states.vx[twallIndex],states.vy[twallIndex],states.x[twallIndex],states.y[twallIndex]);
    }
    t = tnext;
  }

  c1->cd(1);
  TProfile *prof = hist->ProfileX();
  prof->SetTitle("Event Driven Disk Simulation; x-coordinate;proj. density");
  prof->Draw();



  std::cout << "Press ^c to exit" << std::endl;
  theApp.SetIdleTimer(120,".q");
  theApp.Run();
}

System Initialize(){
  std::array<double,N> x;
  std::array<double,N> y;
  std::array<double,N> vx;
  std::array<double,N> vy;
  TRandom3 *r3 = new TRandom3();
  for(int i = 0;i<N;i++){
    std::array<double,DIM> point = generatePoint(x,y);
    x[i] = point[0];
    y[i] = point[1];
    vx[i] = (r3->Rndm()*2*VABS)-VABS;
    vy[i] = pow(100-vx[i]*vx[i],0.5);  //This is part of IC that all vabs are equal
  }
  /** Put in initial conditions
      1. System is at rest on average
      2. All particles have the same absolute velocities at start
  */
  double vxsum = 0.0;
  double vysum = 0.0;
  for(int i =0;i<N;i++){
    vxsum += vx[i];
    vysum += vy[i];
  }
  for(int i = 0;i<N;i++){
    vx[i] = vx[i] - vxsum/N;
    vy[i] = vy[i] - vysum/N;
  }
  System states;
  states.x = x;
  states.y = y;
  states.vx = vx;
  states.vy = vy;
  return states;
}

std::array<double,DIM> generatePoint(std::array<double,N> x,std::array<double,N> y){
  int condition = 0;
  double xtemp = 0.0;  //If debugging, check that were not just getting 
  double ytemp = 0.0; // disks at 0,0
  TRandom3 *r3 = new TRandom3(); //If debugging, check this (two RNG's might be bad)
  while(condition == 0){
    int legalflag = 0;
    xtemp = r3->Rndm()*(SIZE-2*sigma) + sigma;
    ytemp = r3->Rndm()*(SIZE-2*sigma) + sigma;
    for(int i = 0;i<N;i++){ //could make speed improvements here
      if( (fabs(x[i] - 0.0) < 1E-5)&(fabs(y[i]-0.0) <1E-5)){
	continue; //hasn't been initialized yet
      }
      double distance = pow(pow(xtemp-x[i],2) + pow(ytemp-y[i],2),0.5);
      if(distance < 2*sigma){
	legalflag = 1; //Violates space
      }
    }
    if(legalflag == 1){
      condition = 0;
    } else {
      condition = 1;
    }
  }
  std::array<double,DIM> retval = {xtemp,ytemp};
  return retval;
}


std::vector<double> pairCT(System states){
  double tmin = INF; //be careful that you're setting this to an actual time
  double tIndex1 = 0.0;
  double tIndex2 = 0.0;
  for(int i = 1;i<N;i++){  //We don't want to double count or call the function
                           //on itself, i.e. when will it collide with itself
    for(int j = 0;j<i;j++){
      std::array<double,DIM> pairx = {states.x[i],states.x[j]};
      std::array<double,DIM> pairy = {states.y[i],states.y[j]};
      std::array<double,DIM> pairvx = {states.vx[i],states.vx[j]};
      std::array<double,DIM> pairvy = {states.vy[i],states.vy[j]};
      double ti = pairTime(pairx,pairy,pairvx,pairvy);
      if(ti < tmin){
	tmin = ti;
	tIndex1 = j;
	tIndex2 = i;
      }
    }
  }
  std::vector<double> retval;
  retval.push_back(tmin);
  retval.push_back(tIndex1);
  retval.push_back(tIndex2);
  return retval;
}

System pairC(System states,int index1,int index2){
  double delx = states.x[index2]-states.x[index1];
  double dely = states.y[index2]-states.y[index1];
  double abs_x = pow(delx*delx + dely*dely,2);   
  double ePerpx = delx/abs_x;
  double ePerpy = dely/abs_x;
  double delvx = states.vx[index2] - states.vx[index1];
  double delvy = states.vy[index2] - states.vx[index1];
  double scal = delvx*ePerpx + delvy*ePerpy;
  states.vx[index1] += ePerpx*scal;
  states.vy[index1] += ePerpy*scal;
  states.vx[index2] -= ePerpx*scal;
  states.vy[index2] -= ePerpy*scal;
  return states;
}

System wallC(System states,int index,int direction){
  if(direction == 0){
    states.vx[index] = -1*states.vx[index];
  } else {
    states.vy[index] = -1*states.vy[index];
  }
  return states;
}

double pairTime(std::array<double,DIM> pairx,std::array<double,DIM> pairy,
		std::array<double,DIM> pairvx,std::array<double,DIM> pairvy){
  double delx = pairx[1]-pairx[0];  //difference in distance of center in x
  double dely = pairy[1]-pairy[0]; // in y
  double del_SQ = pow(delx,2)+pow(dely,2);  //absolute magnitude in del
  double delvx = pairvx[1] - pairvx[0]; //difference in velocity in x
  double delvy = pairvy[1] - pairvy[0]; // in y 
  double delv_SQ = pow(delvx,2) + pow(delvy,2);  //absolute magnitude of delv
  double scal = delvx*delx + delvy*dely;
  double Upsilon = pow(scal,2) - delv_SQ *(del_SQ - 4.0*pow(sigma,2));
  double deltaT = 0.0;
  if(Upsilon > 0.0 && scal < 0.0){
    deltaT = -1.0*(scal+pow(Upsilon,0.5))/delv_SQ;
  } else {
    deltaT = INF;
  }
  return deltaT;
}


std::vector<double> wallCT(System states){
  double tmin = INF; //be careful that you're setting this to an actual time
  double tIndex = 0.0;
  double direction = 0.0;
  for(int i = 0;i<N;i++){
    double tix = wallTime(states.x[i],states.vx[i]);
    double tiy = wallTime(states.y[i],states.vy[i]);
    double tiMin = std::min(tix,tiy); //smallest collision time for individual
    if(tiMin < tmin){
      tmin = tiMin;
      tIndex = i;
      if(tix < tiy){
	direction = 0.0;
      } else {
	direction = 1.0;
      }
    }
  }
  std::vector<double> retval;
  retval.push_back(tmin);
  retval.push_back(tIndex);
  retval.push_back(direction);
  return retval;
}


double wallTime(double x, double vx){
  double tix = 0.0;
  if(vx > 0.0){
    tix = (SIZE-sigma-x)/vx;
  } else if(vx <0.0){
    tix = (x-sigma)/fabs(vx);
  } else {
    tix = INF;
  }
  return tix;
}


System update(System states,double t){
  for(int i = 0;i<N;i++){
    states.x[i] = states.vx[i]*t + states.x[i];
    states.y[i] = states.vy[i]*t + states.y[i];
  }
  return states;
}
