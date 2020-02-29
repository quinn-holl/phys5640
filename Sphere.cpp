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

const int N = 4; //Number of spheres in simulation
const int TIME = 5000; //How long system will run in # of collisions
const int SUBTIME = 5; //number of updates between collisions
const double sigma = 0.1; //Radius of hard spheres
const double INF = 1000000; //Definition of Infinity here to compare times
const int DIM = 2; //we're only doing 2d Here

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
} states;  //initializes our data structure that we'll use for program

double pairCT(System states);  //outputs time till next pair collision
System pairC(System states);  //updates system after a pair collision
double wallCT(System states);  //outputs time till next wall collision
System wallC(System states);   //updates system after a wall collision
System update(System states,double t); //propogates system through time t
double wallTime(double x, double vx); //one-coordinate calculation
double pairTime(std::array<double,DIM> pairx,std::array<double,DIM> pairy,
		std::array<double,DIM> pairvx,std::array<double,DIM> pairvy);




int main(int argc, char **argv){
  TApplication theApp("App",&argc,argv);
  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  TCanvas *c1 = new TCanvas("c1","Canvas",dw,dh);
  double t = 0; //initialize time variable
  //Need to initialize N hard disks
  states = {};

  for(int i = 0;i++;i<TIME){  //this is in "collision time"
    double tpair = pairCT(states);
    double twall = wallCT(states);
    double tnext = std::min(tpair,twall);
    if(tnext == INF){
      printf("Something's Gone Wrong, Using INF for time step\n");
    }
    for(int m = 1;m<= SUBTIME; m++){
      double deltaT = t +((double)i/m)*(tnext-t);
      states = update(states,deltaT);
      //probably want to sample here for histogram points
    }
    if(twall > tpair){ //updates states to reflect the collision at tnext
      states = pairCT(states);
    } else {
      states = wallCT(states);
    }
    t = tnext;
  }




  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(120,".q");
  theApp.Run();
}


double pairCT(System states){
  double tmin = INF; //be careful that you're setting this to an actual time
  for(int i = 1;i<N;i++){  //We don't want to double count or call the function
                           //on itself, i.e. when will it collide with itself
    for(int j = 0;j<i;j++){
      std::array<double,DIM> pairx = {states[i].x,states[j].x};
      std::array<double,DIM> pairy = {states[i].y,states[j].y};
      std::array<double,DIM> pairvx = {states[i].vx,states[j].vx};
      std::array<double,DIM> pairvy = {states[i].vy,states[j].vy};
      double ti = pairTime(pairx,pairy,pairvx,pairvy);
      tmin = std::min(ti,tmin);
    }
  }
  return tmin;
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


double wallCT(System states){
  double tmin = INF; //be careful that you're setting this to an actual time
  for(int i = 0;i<N;i++){
    double tix = wallTime(states[i].x,states[i].vx);
    double tiy = wallTime(states[i].y,states[i].vy);
    double tiMin = std::min(tix,tiy); //smallest collision time for individual
    tmin = std::min(tiMin,tmin);
  }
  return tmin;
}


double wallTime(double x, double vx){
  double tix = 0.0;
  if(vx > 0.0){
    tix = (1.0-sigma-x)/vx;
  } else if(vx <0.0){
    tix = (x-sigma)/fabs(vx);
  } else {
    tix = INF;
  }
  return tix;
}

System update(System states,double t){
  for(int i = 0;i<N;i++){
    states[i].x = states[i].vx*t;
    states[i].y = states[i].vy*t;
  }
}
