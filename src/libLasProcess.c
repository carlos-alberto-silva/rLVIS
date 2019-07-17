/*########################*/
/*# Functions to process #*/
/*# waveform lidar data  #*/
/*########################*/

/*#######################################*/
/*# Copyright 2006-2016, Steven Hancock #*/
/*# The program is distributed under    #*/
/*# the terms of the GNU General Public #*/
/*# License.    svenhancock@gmail.com   #*/
/*#######################################*/


/*########################################################################*/
/*# This file is part of libCLidar.                                      #*/
/*#                                                                      #*/
/*# libCLidar is free software: you can redistribute it and/or modify    #*/
/*# it under the terms of the GNU General Public License as published by #*/
/*# the Free Software Foundation, either version 3 of the License, or    #*/
/*#  (at your option) any later version.                                 #*/
/*#                                                                      #*/
/*# libCLidar is distributed in the hope that it will be useful,         #*/
/*# but WITHOUT ANY WARRANTY; without even the implied warranty of       #*/
/*#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #*/
/*#   GNU General Public License for more details.                       #*/
/*#                                                                      #*/
/*#    You should have received a copy of the GNU General Public License #*/
/*#    along with libClidar.  If not, see <http://www.gnu.org/licenses/>.#*/
/*########################################################################*/



/*##################################################*/
/*deconvolve using Gold's method*/

float *deconvolve(float *data,int nBins,float **pulse,int pBins,float res,int maxIter,double minChange,char meth)
{
  int i=0,numb=0;
  float *decon=NULL;
  double *deconDo=NULL;
  double *pulseDo=NULL;
  double *dataDo=NULL;
  double *resamplePulse(int,float **,float,int);
  double *goldMeth(double *,double *,int,int,double);
  double *richLucy(double *,double *,int,int,double);
  double energy=0,newEn=0;   /*to balance energies*/


  /*arrays have to be of base 2 length*/
  numb=pow(2.0,(float)((int)(log((double)nBins)/log(2.0)+0.5)+1));
  dataDo=dalloc(numb,"dataDo",0);

  energy=0.0;
  for(i=0;i<nBins;i++){
    dataDo[i]=(double)data[i];
    energy+=dataDo[i];
  }

  /*pulse needs resampling to the correct resolution*/
  pulseDo=resamplePulse(numb,pulse,res,pBins);

  /*call deconvolution method*/
  if(meth==0)     deconDo=goldMeth(dataDo,pulseDo,numb,maxIter,minChange);
  else if(meth==1)deconDo=richLucy(dataDo,pulseDo,numb,maxIter,minChange);
  else{
    fprintf(stderr,"Deconvolution method not defined\n");
    exit(1);
  }
  TIDY(pulseDo);  /*tidy as we go along*/
  TIDY(dataDo);

  /*transfer data*/
  decon=falloc((uint64_t)nBins,"",0);
  newEn=0;
  for(i=0;i<nBins;i++)newEn+=deconDo[i];
  for(i=0;i<nBins;i++)decon[i]=(float)(deconDo[i]*energy/newEn);

  TIDY(deconDo);
  return(decon);
}
/*deconvolve*/


/*##################################################*/
/*Gold's method*/

double *goldMeth(double *data,double *pulse,int numb,int maxIter,double minChange)
{
  int i=0,j=0;
  int real=0,imag=0;
  double *o=NULL,*work=NULL;
  double *smooth=NULL,*denom=NULL;
  double tot=0,changeSq=0,minChangeSq=0;
  double new=0;
  int gsl_fft_complex_radix2_forward(gsl_complex_packed_array,size_t,size_t);
  int gsl_fft_complex_radix2_backward(gsl_complex_packed_array, size_t,size_t);

  minChangeSq=minChange*minChange;

  /*perform FFT*/
  gsl_fft_complex_radix2_forward((gsl_complex_packed_array)pulse,1,numb);

  /*real arrays*/
  o=dalloc(numb,"o",0);
  denom=dalloc(numb,"denominator",0);
  /*complex arrays*/
  work=dalloc(2*numb,"workspace",0);
  smooth=dalloc(2*numb,"workspace",0);

  /*initial guess*/
  for(i=0;i<numb;i++)o[i]=data[i];

  /*iterate over Gold's method*/
  do{
    /*fourier current estimate*/
    for(i=0;i<numb;i++){
      work[2*i]=o[i];
      work[2*i+1]=0.0;
    }
    gsl_fft_complex_radix2_forward((gsl_complex_packed_array)work,1,numb);

    /*blur with pulse*/
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      smooth[real]=pulse[real]*work[real]-pulse[imag]*work[real];
      smooth[imag]=pulse[imag]*work[real]+pulse[real]*work[imag];
    }
    gsl_fft_complex_radix2_backward((gsl_complex_packed_array)smooth,1,numb);

    /*reblur deoniminator*/
    tot=0.0;
    changeSq=0.0;
    for(i=0;i<numb;i++){
      real=i*2;
      imag=i*2+1;
      denom[i]=sqrt(smooth[real]*smooth[real]+smooth[imag]*smooth[imag]);
      if(denom[i]>0.0){
        new=o[i]*data[i]/denom[i];
        tot+=new;
      }else new=0.0;
      changeSq+=(o[i]-new)*(o[i]-new);
      o[i]=new;
    }
    changeSq/=(double)numb;
    /*fprintf(stdout,"Iter %d change %.20f\n",j,changeSq);*/
    if((minChange>=0.0)&&(changeSq<=minChangeSq))break;
    j++;
  }while(j<maxIter);

  TIDY(work);
  TIDY(smooth);
  TIDY(denom);

  return(o);
}/*goldMeth*/




/*##################################################*/
/*process floating point waveform to denoise and deconvolve*/

float *processFloWave(float *wave,int waveLen,denPar *decon,float gbic)
{
  int i=0;
  float *temp=NULL;
  float *mediated=NULL;
  float *preSmoothed=NULL;
  float *mSmoothed=NULL;
  float *denoised=NULL;
  float *smoothed=NULL;
  float *processed=NULL;
  float *denoise(float,float,int,int,float *,float,char);
  float *medianFloat(float *,int,int);
  float *deconvolve(float *,int,float **,int,float,int,double,char);
  float thisTail=0;        /*tail threshold to use here*/
  float *gaussWave=NULL;
  float *fitGaussians(float *,int,denPar *);
  float *CofGhard(float *,uint32_t);
  float *hardHitWave(denPar *,int);
  float *sampled=NULL;
  float *digitise(float *,int,char,float);
  float *correctDrift(float *,int,int,denPar *);
  void medNoiseStats(float *,uint32_t,float *,float *,float *,float,float,char);
  void meanNoiseStats(float *,uint32_t,float *,float *,float *,float,float,int);
  char checkHardEnergy(int *,float *,float);
  char testHard(denPar *,float *,int,float);
  char hardTarget=0;


/*smooth before denoising*/
if(decon->psWidth>0.0){   /*Gaussian smoothing*/
preSmoothed=smooth(decon->psWidth,waveLen,wave,decon->res);
}else if(decon->preMatchF){  /*matched filter*/
preSmoothed=matchedFilter(wave,waveLen,decon,decon->res);
}else{
  preSmoothed=falloc((uint64_t)waveLen,"",0);
  for(i=0;i<waveLen;i++)preSmoothed[i]=wave[i];
}

/*determine noise statistics*/
if(decon->varNoise){
  if(decon->medStats){
    /*convert to a set bit rate*/
    sampled=digitise(preSmoothed,waveLen,decon->bitRate,decon->maxDN);
    medNoiseStats(sampled,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,decon->bitRate);
    TIDY(sampled);
  }else meanNoiseStats(preSmoothed,waveLen,&(decon->meanN),&(decon->thresh),&thisTail,decon->tailThresh,decon->threshScale,(int)(decon->statsLen/decon->res));
}else thisTail=decon->tailThresh;
if(thisTail<0)thisTail=decon->thresh;

/*smooth again if needed*/
if(decon->msWidth>0.0){   /*Gaussian smoothing after noise stats*/
mSmoothed=smooth(decon->msWidth,waveLen,preSmoothed,decon->res);
  TIDY(preSmoothed);
}else{
  mSmoothed=preSmoothed;
  preSmoothed=NULL;
}

/*median filter if needed*/
if(decon->medLen>0){
  mediated=medianFloat(mSmoothed,decon->medLen,waveLen);
  TIDY(mSmoothed);
}else{
  mediated=mSmoothed;
  mSmoothed=NULL;
}

/*correct for detector drift if needed*/
temp=correctDrift(mediated,waveLen,(int)(decon->statsLen/decon->res),decon);
TIDY(mediated);

/*remove background noise*/
if((decon->meanN>0.0)||(decon->thresh>0.0)){
  denoised=denoise(decon->meanN,decon->thresh,decon->minWidth,waveLen,temp,thisTail,decon->noiseTrack);
  TIDY(temp);
}else{
  denoised=temp;
  temp=NULL;
}

/*see if it a single return*/
if(decon->matchHard)hardTarget=testHard(decon,denoised,waveLen,decon->res);
else                hardTarget=0;
if(hardTarget){  /*bunch up the energy*/
processed=CofGhard(denoised,waveLen);
  TIDY(denoised);
  return(processed);
}

/*smooth if required. Note that pulse is smoothed in readPulse()*/
if(decon->sWidth>0.0){
  smoothed=smooth(decon->sWidth,waveLen,denoised,decon->res);
  TIDY(denoised);
}else if(decon->posMatchF){  /*matched filter*/
smoothed=matchedFilter(denoised,waveLen,decon,decon->res);
  TIDY(denoised);
}else{
  smoothed=denoised;
  denoised=NULL;
}

/*scale by GBIC*/
if((gbic>0.0)&&(gbic!=1.0))for(i=0;i<waveLen;i++)smoothed[i]/=gbic;


/*Gaussian fitting*/
if(decon->fitGauss||decon->gaussFilt){
  gaussWave=fitGaussians(smoothed,waveLen,decon);
  /*test for hard target*/
  if(decon->gaussFilt){
    if((decon->nGauss==1)&&(decon->gPar[2]<=decon->hardWidth))hardTarget=1;
    else hardTarget=checkHardEnergy(&decon->nGauss,decon->gPar,decon->hardWidth);
  }else hardTarget=0;

  /*copy Gaussian if to be used*/
  if(hardTarget){ /*single hit*/
  TIDY(gaussWave);
    TIDY(smoothed);
    gaussWave=hardHitWave(decon,waveLen);
  }else if(decon->fitGauss){ /*pass on Gaussian waveform*/
  TIDY(smoothed);
  }else{                /*delete fitting and pass original*/
  TIDY(gaussWave);
    gaussWave=smoothed;
    smoothed=NULL;
  }
}else{   /*don't fit Gaussians*/
  gaussWave=smoothed;
  smoothed=NULL;
  hardTarget=0;
}

/*deconvolve if required*/
if((decon->deconMeth>=0)&&(hardTarget==0)){
  processed=deconvolve(gaussWave,waveLen,decon->pulse,decon->pBins,\
                       decon->res,decon->maxIter,decon->deChang,decon->deconMeth);
  TIDY(gaussWave);
}else{
  processed=gaussWave;    /*otherwise just use the denoised array*/
gaussWave=NULL;
}

return(processed);
}/*processFloWave*/|


/*##################################################*/
/*resample pulse to complex array at correct res*/

double *resamplePulse(int numb,float **pulse,float res,int pBins)
{
  int i=0,bin=0,step=0,maxBin=0;
  float max=0,maxRange=0;
  double *pulseDo=NULL,total=0;

  pulseDo=dalloc(2*numb,"pulseDo",0);
  for(i=2*numb-1;i>=0;i--)pulseDo[i]=0.0;

  /*find the peak to set at zero as pulse is aligned by CofG*/
  max=-1000.0;
  for(i=0;i<pBins;i++){
    if(pulse[1][i]>max){
      max=pulse[1][i];
      maxRange=pulse[0][i];
      maxBin=i;
    }
  }/*peak finding*/

  /*find nearest pulse point to bin, start from centre*/
  step=(int)(res/(pulse[0][1]-pulse[0][0]));
  total=0.0;
  for(i=maxBin;i<pBins;i+=step){
    bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre up*/
  for(i=maxBin-1;i>=0;i-=step){  /*then work from the centre backwards*/
  bin=(int)((pulse[0][i]-maxRange)/res);
    if(bin<0)bin+=numb;
    pulseDo[2*bin]=pulse[1][i];
    total+=pulse[1][i];
  }/*from centre back*/

  /*normalise pulse and prevent zeroes*/
  for(i=0;i<numb;i++){
    if(pulseDo[2*i]<0.0)pulseDo[2*i]=0.0;
    else                pulseDo[2*i]/=total;
  }/*normalisation*/

  return(pulseDo);
}/*reamplePulse*/


void setDenoiseDefault(denPar *denoise)
{
  /*denoising*/
  denoise->meanN=12.0;
  denoise->tailThresh=-1.0;
  denoise->thresh=15.0;
  denoise->minWidth=6;
  denoise->sWidth=0.0;
  denoise->msWidth=0.0;
  denoise->psWidth=0.0;
  denoise->medLen=0;
  denoise->varNoise=0;
  denoise->medStats=0;
  denoise->statsLen=30.0;
  denoise->noiseTrack=1;
  denoise->threshScale=1.0;
  denoise->bitRate=8;
  denoise->maxDN=-1.0;
  denoise->preMatchF=0;    /*no matched filter before denoising*/
  denoise->posMatchF=0;    /*no matched filter after denoising*/
  /*detector drift*/
  denoise->corrDrift=0;    /*do not correct for drift*/
  denoise->varDrift=1;     /*if we do correct, use a variable drift factor*/
  denoise->fixedDrift=0.0; /*if we are using a fixed drift factor, use this*/

  /*deconvolution*/
  denoise->deconMeth=-1;     /*do not deconvolve*/
  denoise->pScale=1.0;      /*scale pulse length by*/
  denoise->maxIter=2000;     /*maximum number of iterations*/
  denoise->deChang=pow(10,0-7.0);  /*change between decon iterations to stop*/
  denoise->deconGauss=1;     /*use Gaussian pulse by default*/
  denoise->pSigma=1.0;
  strcpy(denoise->pNamen,"/home/sh563/data/bess/leica_shape/leicaPulse.dat");  /*pulse filename*/
  denoise->pulse=NULL;       /*pulse to deconvolve by*/
  denoise->pBins=0;          /*number of pulse bins*/
  denoise->res=0.15;

  /*Gaussian fitting*/
  denoise->gWidth=1.5;
  denoise->fitGauss=0;       /*do not fit Gaussians*/
  denoise->gaussPulse=0;     /*do not turn pulse to Gaussian*/
  denoise->minGsig=0.00001;  /*minimum Gaussian fitting width*/

  /*Gaussian hard target identification*/
  denoise->gaussFilt=0;    /*switch*/
  denoise->hardWidth=0.0;  /*maxWidth of hard feature*/
  denoise->hardTol=1.0;    /*tolerance to scale width by*/

  /*correlation hard target finding*/
  denoise->matchHard=0;

  return;
}/*setDenoiseDefault*/
