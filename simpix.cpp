#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TASImage.h"
#include "TApplication.h"
#include "TSystem.h"
#include "assert.h"
#include <iostream>



using namespace std;
/*
** THE PIXELS LIST
** NEED TO CHANGE NMAX DEPENDING ON TRIAL RUN
*/

/*
** Method for making random changes in cities list
*/

int success = 1;
double Distancefind(UInt_t code1,UInt_t code2){
  int Red1 = (UInt_t)((code1 >> 16) & 0xFF);  
  int Red2 = (UInt_t)((code2 >> 16) &0xFF); 
  int Green1 = (UInt_t)((code1 >> 8) & 0xFF);
  int Green2 = (UInt_t)((code2 >> 8) &0xFF);
  int Blue1 = (UInt_t)((code1) & 0xFF);
  int Blue2 =  (UInt_t)((code2) &0xFF);
  double retval = pow((Red1-Red2),2) + pow((Green1-Green2),2) + pow((Blue1-Blue2),2);
  return sqrt(retval);
}

double energy(UInt_t *arr,UInt_t *src,int len);


UInt_t* new_order(UInt_t *arr,int len){
  int size = len; ////*****/////****???
  int rand_1 = (int)(drand48()*(len));
  int rand_2 = (int)((drand48()*(len)));
  UInt_t temp = arr[rand_1];
  arr[rand_1] = arr[rand_2];
  arr[rand_2] = temp;
  return arr;
}

UInt_t* update_order(double beta,UInt_t *arr,int len,UInt_t *src){
  //UInt_t *temp = new_order(arr,len);
  int rand_1 = (int)(drand48()*(len));
  int rand_2 = (int)((drand48()*(len)));
  UInt_t temp = arr[rand_1];
  arr[rand_1] = arr[rand_2];
  arr[rand_2] = temp;
  double DeltaE = 0.0;
  DeltaE = Distancefind(arr[rand_2],src[rand_1])+Distancefind(arr[rand_1],src[rand_2])- Distancefind(arr[rand_1],src[rand_1])-Distancefind(arr[rand_2],src[rand_2]);
  double rand = drand48();
  //printf("%lf\n",DeltaE);
  if(DeltaE <=0 || rand < exp(-beta*DeltaE)){
      UInt_t temp = arr[rand_1];
      arr[rand_1] = arr[rand_2];
      arr[rand_2] = temp;
      success++;
  }
  return arr;
}		    






double energy(UInt_t *arr,UInt_t *src,int len){
  double sum = 0.0;
  for(int i = 0;i<len;i++){
    double iE = Distancefind(arr[i],src[i]);
    sum+= iE;
  }

  return sum;
}


int main(int argc, char *argv[]){
  
  int itherm,ntherm = 9000000;
  double beta,Tmax,T;
  if (argc<3) {
    cout << "Usage: simapix_start image1 image2 <output=out.png>" << endl;
    return 0; 
  }
  TString fsrc=argv[1];
  TString ftgt=argv[2];
  TString fout;
  argc>3 ? fout = argv[3] : fout="out.png";
  cout << "Reading images: source= " << fsrc << " target= " << ftgt << endl;
  cout << "Output= " << fout << endl;
 
  // TApplication theApp("App", &argc, argv);
 
  // create image objects
  TASImage *src = new TASImage(fsrc.Data());
  TASImage *tgt = new TASImage(ftgt.Data());
  TASImage *out = new TASImage(*src); // start with copy of source

  cout << src->GetWidth() << " Src   " << src->GetHeight() << endl;
  cout << tgt->GetWidth() << " TGT   " << tgt->GetHeight() << endl;
 
  // Test image geometry, exit if they are not the same dimensions
  assert ( src->GetWidth() == tgt->GetWidth() && src->GetHeight() == tgt->GetHeight() );
  cout << "Pixel Geometry: " << src->GetWidth() << " x " << src->GetHeight() << endl;
  Long_t numPix=src->GetWidth()*src->GetHeight();
 
  // *** The work happens here
  // access the pixels for the output image 
  // each pixel is a 32-bit word, 1 byte each for (alpha,red,green,blue)
  // don't touch alpha (bits 31:28)
  UInt_t *outPix = out->GetArgbArray();  
  UInt_t *tgtPix = tgt->GetArgbArray();
 
  /*
  // examples of pixel manipulations 
  for (int i=0;i< numPix; i++){
    //  outPix[i]&=0xff00ffff;  // turn off red
    outPix[i]&=0xffff00ff;  // turn off green
    cout << Green(outPix[i]) << endl;
    cout << Red(outPix[i]) << endl;
    //  outPix[i]&=0xffffff00;  // turn off blue
    cout  <<  outPix[i]<<endl;  // print pixel values in hex
  }
  */
  
  srand48((long)time(0));
  T = 8000000;

  
  
  for(int i = 0;i<ntherm;i++){
    beta = 1/T;
    outPix = update_order(beta,outPix,numPix,tgtPix);
  }
  

  while(success>0){
    cout << T << endl;
    success=0;
    beta = 1/T;
    for(itherm=0;itherm<ntherm;itherm++){
      outPix = update_order(beta,outPix,numPix,tgtPix);
    }
    T = T*0.9;
    if(T < 1E-6){
      break;
    }
  }
  
  

  TCanvas *c1 = new TCanvas("c1", "images", 1280, 1024);  //1280x1014
  c1->Divide(2,2);
 
  c1->cd(1);
  c1->Draw();
  src->Draw("X");
  c1->cd(2);
  tgt->Draw("X");
  c1->cd(3);
  out->Draw("X");
  c1->Print("collage.png");
 
  // save the new image
  out->WriteImage(fout.Data());
 
  // cout << "Press ^c to exit" << endl;
  //theApp.Run();

  
  
  return 0;
}
