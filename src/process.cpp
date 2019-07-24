#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "hdf5.h"
extern "C" {
#include "tools/tools.h"
#include "libclidar/libLasRead.h"
#include "libclidar/libLasProcess.h"
#include "libclidar/libLidarHDF.h"
#include "libclidar/libOctree.h"
#include "gedisimulator/gediIO.h"
#include "gedisimulator/gediNoise.h"
#include "gedisimulator/gediMetric.h"

/*#define USEPHOTON*/

  #ifdef USEPHOTON
  #include "photonCount.h"
  #endif

}

#include <Rcpp.h>

  /*##############################*/
  /*# Generates metrics from     #*/
  /*# simulated GEDI waveforms   #*/
  /*# 2015 svenhancock@gmail.com #*/
  /*##############################*/

  /*#######################################*/
  /*# Copyright 2015-2016, Steven Hancock #*/
  /*# The program is distributed under    #*/
  /*# the terms of the GNU General Public #*/
  /*# License.    svenhancock@gmail.com   #*/
  /*#######################################*/


  /*########################################################################*/
  /*# This file is part of the NASA GEDI simulator, gediRat.               #*/
  /*#                                                                      #*/
  /*# gediRat is free software: you can redistribute it and/or modify      #*/
  /*# it under the terms of the GNU General Public License as published by #*/
  /*# the Free Software Foundation, either version 3 of the License, or    #*/
  /*#  (at your option) any later version.                                 #*/
  /*#                                                                      #*/
  /*# gediRat is distributed in the hope that it will be useful,           #*/
  /*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
  /*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
  /*#   GNU General Public License for more details.                       #*/
  /*#                                                                      #*/
  /*#    You should have received a copy of the GNU General Public License #*/
  /*#    along with gediRat.  If not, see <http://www.gnu.org/licenses/>.  #*/
  /*########################################################################*/

  /*element reflectance*/
float rhoG;
float rhoC;

control* makeControl(const char* input, const char* output){
  int i=0;
  control *dimage=NULL;

  /*allocate structures*/
  if(!(dimage=(control *)calloc(1,sizeof(control)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.den=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }
  if(!(dimage->gediIO.gFit=(denPar *)calloc(1,sizeof(denPar)))){
    fprintf(stderr,"error control allocation.\n");
    exit(1);
  }


  /*defaults*/
  /*input/output*/
  dimage->gediIO.nFiles=1;
  dimage->gediIO.inList = chChalloc(dimage->gediIO.nFiles,"inList",0);
  dimage->gediIO.inList[0] = challoc(200,"inList",0);
  strcpy(&(dimage->gediIO.inList[0][0]),"/Users/stevenhancock/data/gedi/analysis/side_lobe/laselva/contLaselva.0.12.wave");
  strcpy(dimage->outRoot,"teastMetric");
  dimage->maxGauss=20;
  dimage->opooMet=NULL;
  dimage->opooGauss=NULL;
  dimage->hdfGedi=NULL;
  /*scan settings*/
  dimage->gediIO.pSigma=0.764331; /*pulse length*/
  dimage->gediIO.fSigma=5.5;      /*footprint width*/
  dimage->gediIO.res=0.15;

  /*switches*/
  dimage->writeFit=0;
  dimage->gediIO.ground=0;
  dimage->gediIO.useInt=0;
  dimage->gediIO.useCount=1;
  dimage->gediIO.useFrac=0;
  dimage->rhRes=5.0;
  dimage->laiRes=10.0;
  dimage->maxLAIh=30.0;
  dimage->bayesGround=0;
  dimage->noise.missGround=0;
  dimage->noise.linkNoise=0;
  dimage->noise.driftFact=0.0;
  dimage->gediIO.linkPsig=0.764331; /*pulse length*/
  dimage->gediIO.linkFsig=5.5;      /*footprint width*/
  dimage->noise.trueSig=5.0;
  dimage->noise.deSig=0.0; //0.1; //4.0*0.15/2.355;
  dimage->noise.bitRate=12;
  dimage->noise.maxDN=4096.0; //1.0/(dimage->pSigma*sqrt(2.0*M_PI));
  dimage->noise.minGap=0.0;
  dimage->noRHgauss=0;    /*do find RH metrics by Gaussian fitting*/
  dimage->renoiseWave=0;  /*do not denoise "truth"*/
  dimage->noise.newPsig=-1.0;   /*leave blank*/
  dimage->gediIO.dontTrustGround=0;  /*do trust ground in waveforms, if there*/
  dimage->readBinLVIS=0;      /*read ASCII rather than binary LVIS*/
  dimage->readHDFlvis=0;      /*read ASCII rather than HDF5 LVIS*/
  dimage->readHDFgedi=0;      /*read ASCII rather than HDF5 GEDI*/
  dimage->gediIO.readPsigma=1;       /*read pSigma from file*/
  dimage->coord2dp=1;         /*round up coords in output*/
  dimage->useBounds=0;        /*process all data provided*/
  dimage->writeGauss=0;       /*do not write Gaussian parameters*/

  /*set default denoising parameters*/
  setDenoiseDefault(dimage->gediIO.den);
  dimage->gediIO.den->meanN=0.0;  /*we haven't added noise yet*/
  dimage->gediIO.den->thresh=0.00000001;  /*tiny number as no noise yet*/
  dimage->gediIO.den->noiseTrack=0;
  dimage->gediIO.den->minWidth=0;
  dimage->gediIO.den->varNoise=0;
  dimage->gediIO.den->threshScale=1.5;
  dimage->gediIO.den->fitGauss=0;
  dimage->gediIO.den->psWidth=0.0;

  /*set default Gaussian fitting  parameters*/
  setDenoiseDefault(dimage->gediIO.gFit);
  dimage->gediIO.gFit->meanN=0.0;    /*no denoising here*/
  dimage->gediIO.gFit->thresh=0.000000005;    /*no denoising here*/
  dimage->gediIO.gFit->noiseTrack=0;    /*no denoising here*/
  dimage->gediIO.gFit->minWidth=0;    /*no denoising here*/
  dimage->gediIO.gFit->varNoise=0;    /*no denoising here*/
  dimage->gediIO.gFit->gWidth=1.2;
  dimage->gediIO.gFit->sWidth=0.0;
  dimage->gediIO.gFit->fitGauss=1;
  dimage->gediIO.gFit->minGsig=0.764331;
  /*noise parameters*/
  dimage->noise.meanN=0.0;
  dimage->noise.nSig=0.0;
  dimage->bThresh=0.001;
  dimage->noise.hNoise=0.0;
  dimage->noise.offset=94.0;
  /*projection, not yet used*/
  dimage->gediIO.wEPSG=4326;  /*waveforms*/
  dimage->gediIO.bEPSG=4326;  /*bounds*/
  /*LVIS data*/
  dimage->lvis.data=NULL;
  dimage->hdfLvis=NULL;
  /*LVIS level2 data*/
  dimage->readL2=0;   /*do not read L2*/
  /*photon counting*/
  dimage->ice2=0;             /*GEDI mode, rather than ICESat-2*/
#ifdef USEPHOTON
  dimage->photonCount.designval=2.1;
  dimage->photonCount.prob=NULL;
  dimage->photonCount.pBins=0;
  dimage->photonCount.H=200.0;
  dimage->photonCount.noise_mult=0.1;
  dimage->photonCount.rhoVrhoG=1.0;
  dimage->photonCount.writeHDF=0;  /*write ASCII by default*/
  dimage->photonCount.hdf=NULL;
#endif
  /*others*/
  rhoG=0.4;
  rhoC=0.57;
  dimage->rhoRatio=rhoC/rhoG;
  dimage->gTol=0.0;
  dimage->gediIO.nMessages=200;
  dimage->fhdHistRes=0.001;

  /*read the command line*/
  TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
  dimage->gediIO.inList=NULL;
  dimage->gediIO.nFiles=1;
  dimage->gediIO.inList=chChalloc(dimage->gediIO.nFiles,"input name list",0);
  dimage->gediIO.inList[0]=challoc((uint64_t)strlen(input)+1,"input name list",0);
  strcpy(dimage->gediIO.inList[0],input);
  strcpy(dimage->outRoot,output);


  /*read deconvolution pulse if needed*/
  if(dimage->gediIO.den->preMatchF||dimage->gediIO.den->preMatchF||dimage->gediIO.den->deconMeth>=0)readPulse(dimage->gediIO.den);
  if((!dimage->gediIO.ground)&&(dimage->noise.missGround)){
    fprintf(stderr,"Noise option conflict. Cannot use missGround without ground\n");
    Rcpp::stop("1");
  }

  return(dimage);
}/*readCommands*/


// [[Rcpp::export]]
Rcpp::NumericVector processFloWave2(Rcpp::CharacterVector input, Rcpp::CharacterVector output)
{
  int i=0;
  control *dimage=NULL;
  dataStruct *data=NULL;
  metStruct *metric=NULL;
  std::string input_str = Rcpp::as<std::string>(input[0]);
  std::string output_str = Rcpp::as<std::string>(output[0]);
  float *processed=NULL,*denoised=NULL;;

  /*read command Line*/
    dimage=makeControl(input_str.c_str(), output_str.c_str());

  /*set link noise if needed*/
    dimage->noise.linkSig=setNoiseSigma(dimage->noise.linkM,dimage->noise.linkCov,dimage->gediIO.linkPsig,dimage->gediIO.linkFsig,rhoC,rhoG);

  /*allocate metric array*/
    if(!(metric=(metStruct *)calloc(1,sizeof(metStruct)))){
      fprintf(stderr,"error metric structure allocation.\n");
      Rcpp::stop("1");
    }

  /*loop over files*/
    for(i=0;i<dimage->gediIO.nFiles;i++){
      if((i%dimage->gediIO.nMessages)==0)Rprintf("Wave %d of %d\n",i+1,dimage->gediIO.nFiles);

    /*read waveform*/
    if(dimage->readBinLVIS)     data=readBinaryLVIS(dimage->gediIO.inList[0],&dimage->lvis,i,&dimage->gediIO);
    else if(dimage->readHDFlvis)data=unpackHDFlvis(dimage->gediIO.inList[0],&dimage->hdfLvis,&dimage->gediIO,i);
    else if(dimage->readHDFgedi)data=unpackHDFgedi(dimage->gediIO.inList[0],&dimage->gediIO,&dimage->hdfGedi,i);
    else                        data=readASCIIdata(dimage->gediIO.inList[i],&(dimage->gediIO));
    if(dimage->readL2)setL2ground(data,i,dimage);

    /*check bounds if needed*/
    if(dimage->useBounds)checkWaveformBounds(data,dimage);

    /*is the data usable*/
    if(data->usable){
      /*denoise and change pulse if needed*/
      if(dimage->renoiseWave)modifyTruth(data,&dimage->noise);

      /*determine truths before noising*/
      determineTruth(data,dimage);

      /*add noise if needed*/
      addNoise(data,&dimage->noise,dimage->gediIO.fSigma,dimage->gediIO.pSigma,dimage->gediIO.res,rhoC,rhoG);

      /*process waveform*/
      /*denoise*/
      denoised=processFloWave(data->noised,data->nBins,dimage->gediIO.den,1.0);

      /*are we in GEDI mode?*/
      if(!dimage->ice2){
        /*Gaussian fit*/
        if(dimage->noRHgauss==0)processed=processFloWave(denoised,data->nBins,dimage->gediIO.gFit,1.0);

        /*shift Gaussian centres to align to absolute elevation*/
        alignElevation(data->z[0],data->z[data->nBins-1],dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss);

        /*determine metrics*/
        findMetrics(metric,dimage->gediIO.gFit->gPar,dimage->gediIO.gFit->nGauss,denoised,data->noised,data->nBins,data->z,dimage,data);

        /*write results*/
        if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)writeResults(data,dimage,metric,i,denoised,processed,dimage->gediIO.inList[0]);
        else                                                             writeResults(data,dimage,metric,i,denoised,processed,dimage->gediIO.inList[i]);
      }else{  /*ICESat-2 mode*/
        photonCountCloud(denoised,data,&dimage->photonCount,dimage->outRoot,i,dimage->gediIO.den,&dimage->noise);
      }/*operation mode switch*/
    }/*is the data usable*/


    /*tidy as we go along*/
    TIDY(processed);
    TIDY(denoised);
    if(data){
      TIDY(data->noised);
      if(dimage->readHDFgedi){  /*pointer to array. do not free*/
        data->wave[0]=NULL;
        if(data->ground)data->ground[0]=NULL;
      }
      TTIDY((void **)data->ground,data->nWaveTypes);
      TTIDY((void **)data->wave,data->nWaveTypes);
      TIDY(data->totE);
      TIDY(data->z);
      TIDY(data);
    }
    TIDY(dimage->gediIO.gFit->gPar);
    TIDY(dimage->gediIO.den->gPar);
    dimage->gediIO.den->nGauss=0;
    dimage->gediIO.gFit->nGauss=0;
    TIDY(metric->rhMax);
    TIDY(metric->rhInfl);
    TIDY(metric->rhReal);
    TIDY(metric->rh);
    TIDY(metric->bGr);
    TIDY(metric->tLAI);
    TIDY(metric->gLAI);
    TIDY(metric->hgLAI);
    TIDY(metric->hiLAI);
    TIDY(metric->hmLAI);
    //TIDY(metric->LmomGau);
    //TIDY(metric->LmomRea);
    //TIDY(metric->LmomInf);
    //TIDY(metric->LmomMax);
  }/*file loop*/

  /*TIDY LVIS data if it was read*/
  if(dimage->readBinLVIS)TIDY(dimage->lvis.data);
  if(dimage->readHDFgedi)dimage->hdfGedi=tidyGediHDF(dimage->hdfGedi);


  if(dimage->writeGauss)Rprintf("Written to %s.gauss.txt\n",dimage->outRoot);
  if(!dimage->ice2)Rprintf("Written to %s.metric.txt\n",dimage->outRoot);
  #ifdef USEPHOTON
  else             Rprintf("Written to %s\n",dimage->photonCount.outNamen);
  #endif


  /*tidy up arrays*/
  tidySMoothPulse();
  TIDY(metric);
  if(dimage){
    if(dimage->lvisL2){
      TIDY(dimage->lvisL2->lfid);
      TIDY(dimage->lvisL2->shotN);
      TIDY(dimage->lvisL2->zG);
      TIDY(dimage->lvisL2);
    }
    if(dimage->readBinLVIS||dimage->readHDFlvis||dimage->readHDFgedi)TTIDY((void **)dimage->gediIO.inList,1);
    else                                        TTIDY((void **)dimage->gediIO.inList,dimage->gediIO.nFiles);
    dimage->gediIO.inList=NULL;
    if(dimage->opooMet){
      fclose(dimage->opooMet);
      dimage->opooMet=NULL;
    }
    if(dimage->opooGauss){
      fclose(dimage->opooGauss);
      dimage->opooGauss=NULL;
    }
    #ifdef USEPHOTON
    if(dimage->photonCount.opoo){
      fclose(dimage->photonCount.opoo);
      dimage->photonCount.opoo=NULL;
    }
    TIDY(dimage->photonCount.prob);
    #endif
    if(dimage->gediIO.den){
      TTIDY((void **)dimage->gediIO.den->pulse,2);
      TIDY(dimage->gediIO.den->matchPulse);
      TIDY(dimage->gediIO.den->hardPulse);
      TIDY(dimage->gediIO.den);
    }
    if(dimage->gediIO.gFit){
      TTIDY((void **)dimage->gediIO.gFit->pulse,2);
      TIDY(dimage->gediIO.gFit);
    }
    dimage->hdfLvis=tidyLVISstruct(dimage->hdfLvis);
    TIDY(dimage);
  }
  return(0);
}/*main*/
