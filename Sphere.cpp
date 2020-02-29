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
}; 

//outputs time till next pair coll
std::vector<double> pairCT(System states,int index1,int index2);  
System pairC(System states);  //updates system after a pair collision
std::vector<double>  wallCT(System states,int index,int direction);  //time till next wall c
System wallC(System states);   //updates system after a wall collision
System update(System states,double t); //propogates system through time t
double wallTime(double x, double vx); //one-coordinate calculation
System Initialize();
double pairTime(std::array<double,DIM> pairx,std::array<double,DIM> pairy,
		std::array<double,DIM> pairvx,std::array<double,DIM> pairvy);




int main(int argc, char **argv){
  TApplication theApp("App",&argc,argv);
  UInt_t dh = gClient->GetDisplayHeight()/2;
  UInt_t dw = 1.1*dh;

  TCanvas *c1 = new TCanvas("c1","Canvas",dw,dh);
  double t = 0; //initialize time variable
  System states = Initialize(); //Need to initialize N hard disks
 

  for(int i = 0;i++;i<TIME){  //this is in "collision time"
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
      double deltaT = t +((double)i/m)*(tnext-t);
      states = update(states,deltaT);
      //probably want to sample here for histogram points
    }
    if(twall > tpair){ //updates states to reflect the collision at tnext
      states = pairCT(states,tpairIndex1,tpairIndex2);
    } else {
      states = wallCT(states,twallIndex,direction);
    }
    t = tnext;
  }




  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(120,".q");
  theApp.Run();
}

System Initialize(){
  std::array<double,N> x = {0.15,0.6,0.43,0.87};
  std::array<double,N> y = {0.65,0.14,0.74,0.0.85};
  std::array<double,N> vx = {4.2,-1.5,3.2,10};
  std::array<double,N> vy = {-5.3,2.6,0.9,10.0};

  /**Idea for initializing in future, fill up array by using rand number 
     in appropriate range (0+sigma,1-sigma) for both x and y, and then checking
     that that's a valid placement. Then do random for vx and vy where we do
     random number between say (-15,15) in both directions
  */
}



std::vector<double> pairCT(System states){
  double tmin = INF; //be careful that you're setting this to an actual time
  double tIndex1 = 0.0;
  double tIndex2 = 0.0;
  for(int i = 1;i<N;i++){  //We don't want to double count or call the function
                           //on itself, i.e. when will it collide with itself
    for(int j = 0;j<i;j++){
      std::array<double,DIM> pairx = {states[i].x,states[j].x};
      std::array<double,DIM> pairy = {states[i].y,states[j].y};
      std::array<double,DIM> pairvx = {states[i].vx,states[j].vx};
      std::array<double,DIM> pairvy = {states[i].vy,states[j].vy};
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
  double delx = states[index2].x-states[index1].x;
  double dely = states[index2].y-states[index1].y;
  double abs_x = pow(delx*delx + dely*dely);   
  double ePerpx = delx/abs_x;
  double ePerpy = dely/abs_x;
  double delvx = states[index2].vx - states[index1].vx;
  double delvy = states[index2].vy - states[index2].vx;
  double scal = delvx*ePerpx + delvy*ePerpy;
  states[index1].x += ePerpx*scal;
  states[index1].y += ePerpy*scal;
  states[index2].x -= ePerpx*scal;
  states[index2].y -= ePerpy*scal;
  return states;
}

System wallC(System states,int index,int direction){
  if(direction == 0){
    states[index].vx = -1*states[index].vx;
  } else {
    states[index].vy = -1*states[index].vy;
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
    double tix = wallTime(states[i].x,states[i].vx);
    double tiy = wallTime(states[i].y,states[i].vy);
    double tiMin = std::min(tix,tiy); //smallest collision time for individual
    if(tiMin < tmin){
      tmin = tiMin;
      tIndex = i;
      if(tix < tiy){
	direction = 0.0;
      } else {
	direction = 1.1;
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
