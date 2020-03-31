#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "TRandom2.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <list>

using namespace std;

const int TRIALS = 5;
double bound;
int n;
int n_steps;

void init_array(double L[][2],double step,double sigma){
  double x = 0;
  double y = 0;
  
  for(int i = 0; i<n; i++){
    x += step;
    if(x - sigma < 0) x += step;
    if(y - sigma < 0) y += step;
    if(x + sigma > bound) y += step;
    if(x + sigma > bound) x = step;
     
    
    L[i][0] = x;
    L[i][1] = y;
  }
  
}


void run_box(double L[][2], double sigma, TH1F *h1){
  double sigma_sq = sigma*sigma;
  double delta = sigma*5;
  int i_event = 0;

  //create cell scheme, cellSize determined by ensurance of no overlap w/ given delta
  double cellSize = 2*sigma + delta;   //check that this is right
  int Xcells = int(bound/cellSize); //number of cells along x-axis 
  cellSize = (1.0/Xcells); //need to recast cellSize to fit the integer steps

  /*  Create Array of Lists, for bookkeeping scheme */
  vector<vector<list<int>>> cells; 
  for(int i = 0;i <Xcells; i++){
    cells.push_back(vector<list<int>>(Xcells));
  }
  
  /* Place disks into list */
  for(int i = 0; i<n ; i++){
    double x = L[i][0];
    double y = L[i][1]; 
    int kx = int(x/cellSize);
    int ky = int(y/cellSize);
    if(kx < Xcells && ky < Xcells){
      cells[kx][ky].push_back(i);
    } else {
      printf("Array Out of Bounds in placing initial config\n");
      exit(1); //Need to terminate program
    }
  }
  
  TRandom *r2 = new TRandom2();

  while(i_event< 5000*n_steps){
    int index = r2->Integer(n);
    double a[2] = {L[index][0],L[index][1]};
    /** Need to get which cell the disk is in */
    int k1[2] = {(int)(a[0]/cellSize), (int)(a[1]/cellSize)};
    /** And which cell it's going to */
    double b[2] = {a[0] + r2->Uniform(-delta,delta), a[1] + r2->Uniform(-delta,delta)};
    int k2[2] = {(int)(b[0]/cellSize), (int)(b[1]/cellSize)};

    vector<int> neighbors; //add all of possible neighbors to a list
    for(int i = k1[0] - 1 ; i <= k1[0] + 1; i++){
      for(int j = k1[1] -1 ; j <= k1[1] + 1; j++){
	if( (i>= 0 && i< Xcells) && (j>= 0 && j<Xcells)){ //otherwise outside of box/grid
	  int disks = cells[i][j].size(); //determines number of disks we need to add
	  if(disks!=0){
	    list<int>::iterator it;
	    for(it = cells[i][j].begin();it!= cells[i][j].end();it++){
	      neighbors.push_back(*it);
	    }
	  }
	}
      }
    }

  
    vector<double> dist; //check all of the neighbors for new distance
    for(unsigned int i = 0; i<neighbors.size();i++){
      if(neighbors[i] != index  && neighbors.size() != 0){
	dist.push_back(pow(b[0]-L[neighbors[i]][0],2) + pow(b[1] - L[neighbors[i]][1],2));
      }
    }
    double min_dist = bound;  //just some large number that will be larger than sigma
    if(dist.size()!=0){
      min_dist = *std::min_element(dist.begin(),dist.end());
    }

    bool box_cond = min(b[0],b[1]) < sigma || max(b[0],b[1]) > bound-sigma;
    if(!(box_cond || min_dist< 4.0*sigma_sq)){
      L[index][0] = b[0];
      L[index][1] = b[1];
      //since we're changing positions, need to update list
      cells[k1[0]][k1[1]].remove(index);
      cells[k2[0]][k2[1]].push_back(index);

      if(i_event % 5000 == 0){
	for(int j = 0;j<n;j++){
	  h1->Fill(L[j][0]);
	  if(L[j][0] > 1-sigma || L[j][0] < sigma){
	    printf("Out of bounds\n");
	  }
	}
      }
      i_event++;
    }
  }
}

void run_markov(double L[][2],double sigma, TH1F *h1){

  double sigma_sq = sigma*sigma;
  double delta = sigma*5;
  int i_event = 0;
  
  TRandom2 *tr = new TRandom2();

  while(i_event < 5000*n_steps){
    int index = tr->Integer(n);
    double a[2]  = {L[index][0],L[index][1]};
    double b[2] = {a[0] + tr->Uniform(-delta,delta),a[1] + tr->Uniform(-delta,delta)};

    vector<double> dist;
    for(int j = 0;j < n; j++){
      if(j != index) dist.push_back(pow(b[0]-L[j][0],2) + pow(b[1]-L[j][1],2));
    }
    
    double min_dist = *std::min_element(dist.begin(),dist.end());
    bool box_cond = min(b[0],b[1]) < sigma || max(b[0],b[1]) > bound - sigma;

    if(!(box_cond || min_dist < 4.0*sigma_sq)){
      L[index][0] = b[0];
      L[index][1] = b[1];
      if(i_event%5000 == 0)
	for(int j = 0;j<n;j++) h1->Fill(L[j][0]);  
      i_event++;
    }
  }
}


int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  n_steps = 10000;
  bound = 1.0;
  n = 64;
  double sigma = 0.05;
  
  double L[n][2];

  //// Loop put disks in grid with as large spacing as possible ////
  double rmax = 0.5;
  bool maximized = false;
  while(!maximized){
    init_array(L,rmax,sigma);
    maximized = true;

    if(rmax < 2*sigma){
      cout<<"not possible with radius"<<endl;
      return 0;
    }
    
    for(int i = 0;i<n;i++)
      if(L[i][0] + sigma > bound || L[i][1] + sigma > bound) maximized = false;
    rmax -= 0.01;
  }

  rmax /= 2; ///This is the maximum radius we can have
  
  TH1F *hMC[TRIALS];
  TH1F *hBox[TRIALS];	
  TCanvas *crMC[TRIALS];
  TCanvas *crBox[TRIALS];	
  
/** we decrement because each successive call to run_markov or run_box
      inherits its initial distribution of disks from the previous runs last
      configuration, so if we ascend sigma, we'll get cases where we begin with
      an invalid configuration to begin with 
  */
  for(int i = 4;i>=0;i--){
    sigma = 0.1*(1+i)*rmax;  ////Set Radius
    TString scanvas1,scanvas2, shist1,shist2,stitle1,stitle2,sfile1,sfile2;
    scanvas1.Form("cBox%d",i);
    scanvas2.Form("cMC%d",i);	  
    sfile1.Form("MCBox_%d.pdf",i);
    sfile2.Form("MC_%d.pdf",i);	  	
    crBox[i] = new TCanvas(scanvas1,"Canvas");
    crMC[i] = new TCanvas(scanvas2,"Canvas");	  
	  
    /** Fill Grid Histograms */
    crBox[i]->cd(1);
    shist1.Form("BoxHist%d",i);	  
    stitle1.Form("MC Density Profile w/ Grid #sigma = %lf;X-val;Count",sigma);
    hBox[i] = new TH1F(shist1,stitle1,200,0,1);
    run_box(L,sigma,hBox[i]);
    hBox[i]->Draw();
    crBox[i]->Update();
    crBox[i]->Print(sfile1);	
	  
    /** Fill MC Histograms */
    crMC[i]->cd(1);
    shist2.Form("MCHist%d",i);	  
    stitle2.Form("MC Density Profile #sigma = %lf;X-val;Count",sigma);
    hMC[i] = new TH1F(shist2,stitle2,200,0,1);
    run_markov(L,sigma,hMC[i]);
    hMC[i]->Draw();
    crMC[i]->Update();
    crMC[i]->Print(sfile2);	  
  }


  
  
  
  cout << "Press ^c to exit" << endl;
  //if(argc == 0){		
  //theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
  //}
}
