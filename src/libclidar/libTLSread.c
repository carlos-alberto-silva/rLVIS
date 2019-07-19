#include <Rinternals.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libLasRead.h"
#include "libLidVoxel.h"


/*######################*/
/*# A library for      #*/
/*# handling tls files #*/
/*# S Hancock, 2016    #*/
/*######################*/


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


/*##################################################################*/
/*write TLS point cloud from binary data*/

void writeTLSpointFromBin(char *namen,double *bounds,FILE *opoo)
{
  uint32_t i=0,nBeams=0;
  uint64_t buffSize=0;   /*buffer size*/
  uint64_t offset=0;
  float zen=0,az=0;
  float xCent=0,yCent=0,zCent=0;
  double x=0,y=0,z=0;
  double xOff=0,yOff=0,zOff=0;
  uint32_t shotN=0;     /*shot number within this scan*/
  uint8_t nHits=0,j=0;  /*number of hits of this beam*/
  float r=0;            /*range*/
  float refl=0;         /*reflectance*/
  char *buffer=NULL;
  FILE *ipoo=NULL;

  /*open file*/
  if((ipoo=fopen(namen,"rb"))==NULL){
    Rprintf("Error opening input file %s\n",namen);
    exit(1);
  }

  /*skip to 4 bytes from the end*/
  if(fseek(ipoo,(long)-4,SEEK_END)){ 
    Rprintf("fseek error from end\n");
    exit(1);
  }
  if(fread(&nBeams,sizeof(uint32_t),1,ipoo)!=1){
    Rprintf("Error reading number of points\n");
    exit(1);
  }
  buffSize=ftell(ipoo);
  /*24 bytes before this are the offset*/
  if(fseek(ipoo,(long)(-1*(4+3*8)),SEEK_END)){ /*skip to 4 bytes from the end*/
    Rprintf("fseek error from end\n");
    exit(1);
  }
  buffer=challoc(3*sizeof(double),"buffer",0);
  if(fread(&(buffer[0]),sizeof(double),3,ipoo)!=3){
    Rprintf("Error reading 3 offsets\n");
    exit(1);
  }
  memcpy(&xOff,&buffer[0],8);
  memcpy(&yOff,&buffer[8],8);
  memcpy(&zOff,&buffer[16],8);
  TIDY(buffer);


  /*read data into buffer*/
  if(fseek(ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
    Rprintf("fseek error to start\n");
    exit(1);
  }
  buffer=challoc((uint64_t)buffSize,"buffer",0);   /*allocate spave*/
  if(fread(&buffer[0],sizeof(char),buffSize,ipoo)!=buffSize){  /*read beams*/
    Rprintf("Error reading point data for writing\n");
    exit(1);
  }
  /*close file*/
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }


  /*loop along buffer and output*/
  offset=0;
  for(i=0;i<nBeams;i++){
    memcpy(&zen,&buffer[offset],4);
    zen*=M_PI/180.0; /*convert to radians*/
    offset+=4;
    memcpy(&az,&buffer[offset],4);
    az*=M_PI/180.0; /*convert to radians*/
    offset+=4;
    memcpy(&xCent,&buffer[offset],4);
    offset+=4;
    memcpy(&yCent,&buffer[offset],4);
    offset+=4;
    memcpy(&zCent,&buffer[offset],4);
    offset+=4;
    memcpy(&shotN,&buffer[offset],4);
    offset+=4;
    memcpy(&nHits,&buffer[offset],1);
    offset+=1;

    /*loop over hits in this beam*/
    for(j=0;j<nHits;j++){
      memcpy(&r,&buffer[offset],4);
      offset+=4;
      memcpy(&refl,&buffer[offset],4);
      offset+=4;

      x=(double)r*sin((double)zen)*cos((double)az)+(double)xCent+xOff;
      y=(double)r*sin((double)zen)*sin((double)az)+(double)yCent+yOff;
      z=(double)r*cos((double)zen)+(double)zCent+zOff;

      /*check bounds*/
      if((x>=bounds[0])&&(y>=bounds[1])&&(z>=bounds[2])&&(x<=bounds[3])&&(y<=bounds[4])&&(z<=bounds[5])){
        fprintf(opoo,"%.3f %.3f %.3f %f %d %d %f %f %f\n",x,y,z,refl,j,nHits,zen,az,r);
      }
    }/*hit loop*/

  }/*beam loop*/
  TIDY(buffer);
  return;
}/*writeTLSpointFromBin*/


/*##################################################################*/
/*read single TLS scan, all data*/

void readTLSpolarBinary(char *namen,uint32_t place,tlsScan **scan)
{
  uint32_t i=0;
  uint64_t offset=0;
  uint64_t buffSize=0;   /*buffer size*/
  uint64_t buffEnd=0;    /*end of file reading buffer*/
  unsigned char j=0;
  char *buffer=NULL;

  /*file reading buffer size. Hardwired to 500 Mbytes for now*/
  buffSize=500000000;

  /*is this the first call?*/
  if((*scan)==NULL){  /*if so, read size and allocate space*/
    Rprintf("Reading %s ",namen);

    /*allocate structure and open file*/
    if(!((*scan)=(tlsScan *)calloc(1,sizeof(tlsScan)))){
      Rprintf("error scan allocation.\n");
      exit(1);
    }
    if(((*scan)->ipoo=fopen(namen,"rb"))==NULL){
      Rprintf("Error opening input file %s\n",namen);
      exit(1);
    }
  
    /*last 4 bytes state number of beams*/
    if(fseek((*scan)->ipoo,(long)-4,SEEK_END)){ /*skip to 4 bytes from the end*/
      Rprintf("fseek error from end\n");
      exit(1);
    }
    (*scan)->totSize=(uint64_t)ftell((*scan)->ipoo);
    if(fread(&((*scan)->nBeams),sizeof(uint32_t),1,(*scan)->ipoo)!=1){
      Rprintf("Error reading number of points\n");
      exit(1);
    }
    Rprintf("There are %d TLS beams\n",(*scan)->nBeams);
    if(buffSize>(*scan)->totSize)buffSize=(*scan)->totSize;  /*reset buffer if it is too big*/
    /*24 bytes before this are the offset*/
    if(fseek((*scan)->ipoo,(long)(-1*(4+3*8)),SEEK_END)){ /*skip to 4 bytes from the end*/
      Rprintf("fseek error from end\n");
      exit(1);
    }
    buffer=challoc(3*8,"buffer",0);
    if(fread(buffer,8,3,(*scan)->ipoo)!=3){
      Rprintf("Error reading 3 offsets\n");
      exit(1);
    }
    memcpy(&((*scan)->xOff),&buffer[0],8);
    memcpy(&((*scan)->yOff),&buffer[8],8);
    memcpy(&((*scan)->zOff),&buffer[16],8);
    TIDY(buffer);

    /*seek back to start of file and set buffer sizes*/
    if(fseek((*scan)->ipoo,(long)0,SEEK_SET)){ /*rewind to start of file*/
      Rprintf("fseek error to start\n");
      exit(1);
    }
    (*scan)->totRead=0;
    (*scan)->pOffset=0;
    (*scan)->nRead=0;
    (*scan)->maxRead=(uint32_t)(buffSize/(sizeof(tlsBeam)-2*4));  /*if all beams had hits, there would be fewer bytes per beam*/
    /*allocate space for beams*/
    if(!((*scan)->beam=(tlsBeam *)calloc((long)(*scan)->maxRead,sizeof(tlsBeam)))){
      Rprintf("error beam allocation. Allocating %lu\n",buffSize);
      exit(1);
    }
  }else if((place==0)||(((uint64_t)place-(uint64_t)(*scan)->pOffset)<(uint64_t)(*scan)->nRead)){  /*do we need to read anymore?*/
    return;  /*if not, pass straight back to avoid reading data*/
  }

  /*allocate buffer space and read. Fudge to prevent reading off the end of the file*/
  if((buffSize+(*scan)->totRead)>(*scan)->totSize)buffSize=(*scan)->totSize-(*scan)->totRead;  /*adjust if at file end*/
  buffer=challoc(buffSize,"buffer",0);      /*allocate space*/
  if(fread(&buffer[0],sizeof(char),buffSize,(*scan)->ipoo)!=buffSize){  /*read beams*/
    Rprintf("Error reading beam data for buffer of size %lu\n",buffSize);
    exit(1);
  }

  /*free up old space to prevent it being reallocated and wasting RAM*/
  for(i=0;i<(*scan)->maxRead;i++){
    TIDY((*scan)->beam[i].r);
    TIDY((*scan)->beam[i].refl);
    (*scan)->beam[i].nHits=0;
  }

  /*load buffer into structure*/
  i=0;
  offset=0;
  buffEnd=buffSize-(2*4+3*8+4+1); /*the longest buffer for a beam with 20 hits*/
  while((offset<buffEnd)&&(i<(*scan)->maxRead)){
    /*copy over a beam*/
    memcpy(&((*scan)->beam[i].zen),&buffer[offset],4);
    (*scan)->beam[i].zen*=M_PI/180.0; /*convert to radians*/
    offset+=4;
    memcpy(&((*scan)->beam[i].az),&buffer[offset],4);
    (*scan)->beam[i].az*=M_PI/180.0; /*convert to radians*/
    offset+=4;
    memcpy(&(*scan)->beam[i].x,&buffer[offset],4);
    offset+=4;
    memcpy(&(*scan)->beam[i].y,&buffer[offset],4);
    offset+=4;
    memcpy(&(*scan)->beam[i].z,&buffer[offset],4);
    offset+=4;
    memcpy(&((*scan)->beam[i].shotN),&buffer[offset],4);
    offset+=4;
    memcpy(&((*scan)->beam[i].nHits),&buffer[offset],1);
    offset+=1;

    /*check we're not going to run off the end*/
    if((offset+8*(uint64_t)(*scan)->beam[i].nHits)>=(*scan)->totSize){
      offset-=2*4+3*8+4+1;  /*rewind counters*/
      i--;
      break;
    }

    /*copy hit information, multiple per beam*/
    if((*scan)->beam[i].nHits>0){
      (*scan)->beam[i].r=falloc((uint64_t)(*scan)->beam[i].nHits,"range",i);
      (*scan)->beam[i].refl=falloc((uint64_t)(*scan)->beam[i].nHits,"refl",i);
      for(j=0;j<(*scan)->beam[i].nHits;j++){
        memcpy(&((*scan)->beam[i].r[j]),&buffer[offset],4);
        offset+=4;
        memcpy(&((*scan)->beam[i].refl[j]),&buffer[offset],4);
        offset+=4;
      }/*hit loop*/
    }/*hit check*/

    i++;  /*increment beam counter*/
  }/*load into array loop*/

  /*increment position*/
  (*scan)->pOffset+=(*scan)->nRead; /*offset to first point in RAM*/
  (*scan)->nRead=i;               /*number of points read this time*/
  (*scan)->totRead+=offset;       /*file positrion in bytes*/

  /*if this is the last call, close file*/
  if(((*scan)->pOffset+(*scan)->nRead)>=(*scan)->nBeams){
    if((*scan)->ipoo){
      fclose((*scan)->ipoo);
      (*scan)->ipoo=NULL;
    }
  }/*file closing check*/

  /*We will not have read a whole number of beams,*/
  /*as each can have a variable number of hits*/
  if((*scan)->ipoo){
    if(fseek((*scan)->ipoo,(long)((*scan)->totRead),SEEK_SET)){ /*rewind to whole number of beams*/
      Rprintf("fseek error to start\n");
      exit(1);
    }
  }

  TIDY(buffer);
  return;
}/*readTLSpolarBinary*/


/*###################################*/
/*read a PTX file*/

/*###################################*/
/*read a PTX file*/

void readPTXleica(char *namen,uint32_t place,tlsScan **scan)
{
  int j=0;
  uint32_t i=0;
  int64_t xInd=0,yInd=0;
  uint32_t contN=0;
  uint64_t fStart=0;
  double x=0,y=0,z=0;
  float zen=0,az=0,diff=0;
  float lastZen=-400.0,lastAz=-400.0;
  char line[100];
  char temp1[25],temp2[25];
  char temp3[25],temp4[25];
  void translateLeica(double *,double *,double *,float **);
  /*keep a list of done beams, only in here*/
  static uint32_t nCol,nRow;
  static double res;
  static float minZen,minAz;

  /*read header if first time*/
  if((*scan)==NULL){
    res=0.0;
    contN=0;
    /*allocate space*/
    if(!((*scan)=(tlsScan *)calloc(1,sizeof(tlsScan)))){
      Rprintf("error in tls structure allocation.\n");
      exit(1);
    }
    (*scan)->matrix=fFalloc(4,"translation matrix",0);
    for(j=0;j<4;j++)(*scan)->matrix[j]=falloc(4,"translation matrix",j+1);

    /*open file*/
    if(((*scan)->ipoo=fopen(namen,"r"))==NULL){
      Rprintf("Error opening input file %s\n",namen);
      exit(1);
    }

    /*read header and determine file length*/
    i=0;
    minAz=minZen=1000.0;
    while(fgets(line,100,(*scan)->ipoo)!=NULL){
      /*read the header*/
     if(i<10){
        if(i==0)nCol=atoi(line);
        else if(i==1)nRow=atoi(line);
        else{
          if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){  /*translation matrix*/
            (*scan)->matrix[0][j]=atof(temp1);
            (*scan)->matrix[1][j]=atof(temp2);
            (*scan)->matrix[2][j]=atof(temp3);
            (*scan)->matrix[3][j]=atof(temp4);
            j++;
          }else if(sscanf(line,"%s %s %s",temp1,temp2,temp3)==3){
            if(i==2){
              (*scan)->xOff=atof(temp1);
              (*scan)->yOff=atof(temp2);
              (*scan)->zOff=atof(temp3);
              j=0;
            }
          }
        }
      }else if(i==10)fStart=ftell((*scan)->ipoo);
      else{  /*work out resolution*/
        if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
          x=atof(temp1);
          y=atof(temp2);
          z=atof(temp3);
          if((fabs(x)+fabs(y)+fabs(z))>0.0001){
            translateLeica(&x,&y,&z,(*scan)->matrix);
            zen=atan2(sqrt(x*x+y*y),z);
            az=atan2(x,y);
            if((lastZen>=-M_PI)&&(lastZen<M_PI)){
              diff=sqrt(pow(lastZen-zen,2)+pow(lastAz-az,2));
              if(diff<1.0){  /*avoid ends of scan lines*/
                res+=(double)diff;
                contN++;
              }
            }
            lastZen=zen;
            lastAz=az;
            /*keep track of min angles*/
            if(zen<minZen)minZen=zen;
            if(az<minAz)minAz=az;
          }else lastZen=lastAz=-400.0;
        }
      }
      i++;
    }
    res/=(double)contN;
    (*scan)->nBeams=i-(10+1); /*ignore header*/  /*NOTE: I am not sure why it is -11 rather than -10, but 10 causes a seg fault*/
    (*scan)->maxRead=500000000/((uint64_t)sizeof(tlsBeam)+8);
    (*scan)->nRead=0;
    (*scan)->totRead=0;
    (*scan)->pOffset=0;
    Rprintf("Scan contains %u beams\n",(*scan)->nBeams);

    /*allocate space*/
    if(!((*scan)->beam=(tlsBeam *)calloc((long)(*scan)->maxRead,sizeof(tlsBeam)))){
      Rprintf("error beam allocation.\n");
      exit(1);
    }
    /*rewind to start of data blocks*/
    if(fseek((*scan)->ipoo,(long)fStart,SEEK_SET)){
      Rprintf("fseek error to start\n");
      exit(1);
    }
  }else if((place-(*scan)->pOffset)<(*scan)->nRead){  /*if still within block, return*/
    return;
  }/*header reading*/


  /*update offset*/
  (*scan)->pOffset+=(*scan)->nRead;

  /*clear out old memory*/
  for(i=0;i<(*scan)->maxRead;i++){
    TIDY((*scan)->beam[i].refl);
    TIDY((*scan)->beam[i].r);
    (*scan)->beam[i].nHits=0;
  }

  /*read block of data*/
  i=0;
  while(fgets(line,100,(*scan)->ipoo)!=NULL){
    if(sscanf(line,"%s %s %s %s",temp1,temp2,temp3,temp4)==4){
      x=atof(temp1);
      y=atof(temp2);
      z=atof(temp3);

      (*scan)->beam[i].x=0.0;  /*origin is always 0 for a Leica*/
      (*scan)->beam[i].y=0.0;
      (*scan)->beam[i].z=0.0;
      (*scan)->beam[i].shotN=i;

      /*contains hits?*/
      if((fabs(x)+fabs(y)+fabs(z))>0.0){
        /*get angles*/
        translateLeica(&x,&y,&z,(*scan)->matrix);
        (*scan)->beam[i].zen=(float)atan2(sqrt(x*x+y*y),z);
        (*scan)->beam[i].az=(float)atan2(x,y);
        /*mark hit*/
        (*scan)->beam[i].nHits=1;
        (*scan)->beam[i].r=falloc(1,"range",i+1);
        (*scan)->beam[i].refl=falloc(1,"refl",i+1);
        (*scan)->beam[i].r[0]=sqrt(x*x+y*y+z*z);
        (*scan)->beam[i].refl[0]=atof(temp4);

        /*mark which pixels have been done for gap tracking*/
        xInd=(int64_t)((double)((*scan)->beam[i].az-minAz)/res+0.5);
        yInd=(int64_t)((double)((*scan)->beam[i].zen-minZen)/res+0.5);
      }else{  /*otherwise gap. NOTE that this should be diaganol across x*/
        if(yInd>0)yInd--;
        else{
          yInd=(int64_t)nRow;
          if(xInd<(int64_t)nCol)xInd++;
          else                  xInd=0;
        }
        (*scan)->beam[i].zen=(float)(res*(double)yInd+(double)minZen);
        (*scan)->beam[i].az=(float)(res*(double)xInd+(double)minAz);
        (*scan)->beam[i].nHits=0;
        (*scan)->beam[i].shotN=i;
        (*scan)->beam[i].r=NULL;
        (*scan)->beam[i].refl=NULL;
      }

      /*the Leica can have angles greater than 180. Wrap around to avoid breaking functions*/
      if((*scan)->beam[i].az<M_PI)(*scan)->beam[i].az+=2.0*M_PI;
      if((*scan)->beam[i].az>M_PI)(*scan)->beam[i].az-=2.0*M_PI;
      if((*scan)->beam[i].zen<M_PI)(*scan)->beam[i].zen+=2.0*M_PI;
      if((*scan)->beam[i].zen>M_PI)(*scan)->beam[i].zen-=2.0*M_PI;

      i++;
      /*if we have read the whole chunk, break*/
      if(i>=(*scan)->maxRead){
        (*scan)->totRead=ftell((*scan)->ipoo);
        break;
      }
    }
  }/*data reading loop*/
  (*scan)->nRead=i;

  Rprintf("Read %u of %u beams\n",(*scan)->pOffset,(*scan)->nBeams);


  /*if this is the last call, close file*/
  if(((*scan)->pOffset+(*scan)->nRead)>=(*scan)->nBeams){
    if((*scan)->ipoo){
      fclose((*scan)->ipoo);
      (*scan)->ipoo=NULL;
    }
    TTIDY((void **)(*scan)->matrix,4);
    (*scan)->matrix=NULL;
  }/*file closing check*/

  return;
}/*readPTXleica*/


/*###################################*/
/*translate leica coord*/

void translateLeica(double *x,double *y,double *z,float **matrix)
{
  double tempX=0,tempY=0,tempZ=0;

  tempX=*x;
  tempY=*y;
  tempZ=*z;

  *x=tempX*matrix[0][0]+tempY*matrix[1][0]+tempZ*matrix[2][0];
  *y=tempX*matrix[0][1]+tempY*matrix[1][1]+tempZ*matrix[2][1];
  *z=tempX*matrix[0][2]+tempY*matrix[1][2]+tempZ*matrix[2][2];

  return;
}/*translateLeica*/


/*##################################################################*/
/*add up gap fraction for intersected voxels*/

void noteVoxelGaps(int *voxList,int nIntersect,double *rangeList,voxStruct *vox,tlsScan *tempTLS,uint32_t j,float maxR,char useFracGap,lidVoxPar *lidPar,int fInd)
{
  int k=0,n=0;
  float lastHitR=0;
  float appRefl=0,rad=0;
  char doIt=0,hasHit=0;


  /*find the range to the last hit*/
  if(tempTLS->beam[j].nHits>0)lastHitR=(tempTLS->beam[j].r[tempTLS->beam[j].nHits-1]<maxR)?tempTLS->beam[j].r[tempTLS->beam[j].nHits-1]:maxR;
  else                        lastHitR=maxR;

  /*loop over intersected voxels*/
  for(k=0;k<nIntersect;k++){
    /*hits before voxel*/
    if(!useFracGap){  /*simple method. All hit until last return*/
      if(rangeList[k]<=lastHitR){  /*entry point is before last return*/
        vox->hits[fInd][voxList[k]]+=1.0;
        doIt=1;
      }else{                       /*entry point is after last return*/
        vox->miss[fInd][voxList[k]]+=1.0;
        doIt=0;
      }/*hit to voxel check*/
    }else{            /*John's fractional method*/
      Rprintf("John's folly method not implemented yet\n");
      exit(1);
    }/*hits before voxel*/

    /*add up total length of beams passing through*/
    vox->totVol[fInd][voxList[k]]+=rangeList[k+1]-rangeList[k];

    /*hits within voxel*/
    if(doIt){  /*only if the beam has made it this far*/
      /*loop over all hits along beam to see which are within voxel*/
      hasHit=0;
      appRefl=0.0;
      for(n=0;n<tempTLS->beam[j].nHits;n++){/*hit along beam loop*/
        if(!lidPar->correctR)appRefl+=(float)tempTLS->beam[j].refl[n];   /*total reflectance to account for occlusion*/
        else                 appRefl+=(float)tempTLS->beam[j].refl[n]*pow((float)tempTLS->beam[j].r[n],2.0);
        if((tempTLS->beam[j].r[n]>=rangeList[k])&&(tempTLS->beam[j].r[n]<rangeList[k+1])){
          hasHit=1;
          /*count up area of points within voxel*/
          if(lidPar){
            rad=tlsPointSize((double)tempTLS->beam[j].r[n],tempTLS->beam[j].refl[n],lidPar->beamTanDiv,\
                                lidPar->beamRad,lidPar->minRefl,lidPar->maxRefl,lidPar->appRefl,1.0);
            vox->sumRsq[fInd][voxList[k]]+=rad*rad;
          }/*count up area of points within voxel*/
        }else if(tempTLS->beam[j].r[n]>=rangeList[k+1])break;  /*left the voxel*/
      }/*hit along beam loop*/
      /*count up number of hits and misses within voxel*/
      if(hasHit){
        vox->inHit[fInd][voxList[k]]+=1.0;
        /*count up volume sampled*/
        if(n<(tempTLS->beam[j].nHits-1)){
          vox->sampVol[fInd][voxList[k]]+=rangeList[k+1]-rangeList[k];  /*not last return*/
        }else{
          vox->sampVol[fInd][voxList[k]]+=tempTLS->beam[j].r[tempTLS->beam[j].nHits-1]-rangeList[k];  /*last return*/
        }
        /*normalise mean refl for waveform gap fraction*/
        if(lidPar){
          appRefl/=(float)(lidPar->maxRefl-lidPar->minRefl);
          vox->meanRefl[fInd][voxList[k]]+=appRefl;
        }
      }else{ /*no hits in this voxel*/
        vox->inMiss[fInd][voxList[k]]+=1.0;
        vox->sampVol[fInd][voxList[k]]+=rangeList[k+1]-rangeList[k];
      }/*hit within voxel check*/
      /*keep track of mean zenith*/
      vox->meanZen[fInd][voxList[k]]+=tempTLS->beam[j].zen;
    }/*beam made it to voxel check*/
  }/*voxel intersection loop*/

  return;
}/*noteVoxelGaps*/


/*##################################################################*/
/*save relevant TLS points*/

void saveTLSpoints(tlsScan *tempTLS,uint32_t j,voxStruct *vox,tlsScan *scan,double xCent,double yCent,double zCent,tlsVoxMap *map,int fInd)
{
  int k=0,vPlace=0;
  int xBin=0,yBin=0,zBin=0;
  int *markInt(int,int *,int);
  uint32_t *markUint32(int,uint32_t *,uint32_t);
  double x=0,y=0,z=0;

  /*loop over hits in this beam*/
  for(k=0;k<tempTLS->beam[j].nHits;k++){
    x=xCent+tempTLS->beam[j].r[k]*sin((double)tempTLS->beam[j].az)*sin((double)tempTLS->beam[j].zen);
    y=yCent+tempTLS->beam[j].r[k]*cos((double)tempTLS->beam[j].az)*sin((double)tempTLS->beam[j].zen);
    z=zCent+tempTLS->beam[j].r[k]*cos((double)tempTLS->beam[j].zen);

    /*check bounds and copy point if within*/
    if((x>=vox->bounds[0])&&(y>=vox->bounds[1])&&(z>=vox->bounds[2])&&\
       (x<=vox->bounds[3])&&(y<=vox->bounds[4])&&(z<=vox->bounds[5])){
      /*voxel space coordinates*/
      xBin=(int)((x-vox->bounds[0])/vox->res[0]);
      yBin=(int)((y-vox->bounds[1])/vox->res[1]);
      zBin=(int)((z-vox->bounds[2])/vox->res[2]);
      vPlace=zBin*vox->nX*vox->nY+yBin*vox->nX+xBin;

      /*are we within voxel space, to avoid rounding errors*/
      if((xBin<0)||(xBin>=vox->nX)||(yBin<0)||(yBin>=vox->nY)||(zBin<0)||(zBin>=vox->nZ)){
        continue;
      }

      /*mark TLS points*/
      scan->point[scan->nPoints].x=(float)(x-scan->xOff);  /*subtract offset to save disk space*/
      scan->point[scan->nPoints].y=(float)(y-scan->yOff);  /*subtract offset to save disk space*/
      scan->point[scan->nPoints].z=(float)(z-scan->zOff);  /*subtract offset to save disk space*/
      scan->point[scan->nPoints].r=tempTLS->beam[j].r[k];
      scan->point[scan->nPoints].refl=tempTLS->beam[j].refl[k];
      scan->point[scan->nPoints].hitN=k;
      scan->point[scan->nPoints].nHits=tempTLS->beam[j].nHits;
      /*map to voxels*/
      map->mapFile[vPlace]=markInt(map->nIn[vPlace],&(map->mapFile[vPlace][0]),fInd);
      map->mapPoint[vPlace]=markUint32(map->nIn[vPlace],&(map->mapPoint[vPlace][0]),scan->nPoints);
      map->nIn[vPlace]++;
      scan->nPoints++;
    }
  }/*hit in beam loop*/

  return;
}/*saveTLSpoints*/


/*##################################################################*/
/*read a single TLS scan within a voxel*/

tlsScan *readOneTLS(char *namen,voxStruct *vox,char useFracGap,tlsVoxMap *map,int fInd,lidVoxPar *lidPar)
{
  int k=0;
  int pInd=0;
  int nBuff=0,vPlace=0;
  int *voxList=NULL,nIntersect=0;
  uint32_t j=0,tInd=0;
  float maxR=0;
  double grad[3],*rangeList=NULL;
  double xCent=0,yCent=0,zCent=0;
  tlsScan *scan=NULL,*tempTLS=NULL;
  void noteVoxelGaps(int *,int,double *,voxStruct *,tlsScan *,uint32_t,float,char,lidVoxPar *,int);
  void saveTLSpoints(tlsScan *,uint32_t,voxStruct *,tlsScan *,double,double,double,tlsVoxMap *,int);
  char checkIfPtx(char *);
  char isPtx=0;  /*ptx file flag*/

  /*is this a ptx file?*/
  isPtx=checkIfPtx(namen);

  /*max range of Riegl: OTHERS ARE SHORTER or longer. COULD BE ADJUSTABLE*/
  maxR=300.0;

  /*allocate space*/
  if(!(scan=(tlsScan *)calloc(1,sizeof(tlsScan)))){
    Rprintf("error in tls structure allocation.\n");
    exit(1);
  }

  /*read all data into RAM*/
  if(isPtx==0)readTLSpolarBinary(namen,0,&tempTLS);
  else        readPTXleica(namen,0,&tempTLS);

  /*if we are saving points, allocate a buffer*/
  if(vox->savePts){
    Rprintf("We will be saving points\n");
    nBuff=4*tempTLS->nBeams;
    if(!(scan->point=(tlsPoint *)calloc(nBuff,sizeof(tlsPoint)))){
      Rprintf("error in tls point allocation. Allocating %lu\n",nBuff*sizeof(tlsPoint));
      exit(1);
    }
    if(map->mapFile==NULL){
      map->mapFile=iIalloc(vox->nVox,"voxel file map",0);        /*file per voxel*/
      if(!(map->mapPoint=(uint32_t **)calloc(vox->nVox,sizeof(uint32_t *)))){
        Rprintf("error in voxel point map allocation.\n");
        exit(1);
      }
      map->nIn=ialloc(vox->nVox,"voxel map number",0);        /*file per voxel*/
    }
  }
  scan->nPoints=scan->nBeams=0;
  scan->beam=NULL;
  scan->xOff=tempTLS->xOff;
  scan->yOff=tempTLS->yOff;
  scan->zOff=tempTLS->zOff;

  /*is the scan within maxR of the bounds?*/
  if(((vox->bounds[0]-tempTLS->xOff)<=maxR)&&((vox->bounds[1]-tempTLS->yOff)<=maxR)&&\
     ((vox->bounds[2]-tempTLS->zOff)<=maxR)&&((vox->bounds[3]-tempTLS->xOff)>=(-1.0*maxR))&&\
     ((vox->bounds[4]-tempTLS->yOff)>=(-1.0*maxR))&&((vox->bounds[5]-tempTLS->zOff)>=(-1.0*maxR))){

    /*loop over beams in scan*/
    for(j=0;j<tempTLS->nBeams;j++){
      /*update where we are in the file if needed*/
      if(isPtx==0)readTLSpolarBinary(namen,j,&tempTLS);
      else        readPTXleica(namen,j,&tempTLS);
      tInd=j-tempTLS->pOffset;   /*update index to account for buffered memory*/

      /*avoid tilt mount if needed*/
      if(fabs(tempTLS->beam[tInd].zen)>=vox->maxZen)continue;  /*skip if zenith too high*/

      /*apply offset to centre*/
      xCent=(double)tempTLS->beam[tInd].x+tempTLS->xOff;
      yCent=(double)tempTLS->beam[tInd].y+tempTLS->yOff;
      zCent=(double)tempTLS->beam[tInd].z+tempTLS->zOff;

      /*check angles are sensible*/
      /*if((tempTLS->beam[tInd].zen<-400.0)||(tempTLS->beam[tInd].zen>400.0)||\
              (tempTLS->beam[tInd].az<-400.0)||(tempTLS->beam[tInd].az>400.0))continue;*/

      /*find intersecting voxels*/
      grad[0]=tempTLS->beam[tInd].zen;
      grad[1]=tempTLS->beam[tInd].az;
      grad[2]=-99999.0;
      voxList=findVoxels(&(grad[0]),xCent,yCent,zCent,vox->bounds,\
                  &(vox->res[0]),&nIntersect,vox->nX,vox->nY,vox->nZ,&rangeList);
      if(nIntersect==0)continue;   /*if no voxels intersected*/

      /*add up gap fraction for intersected voxels*/
      noteVoxelGaps(voxList,nIntersect,rangeList,vox,tempTLS,tInd,maxR,useFracGap,lidPar,fInd);

      /*record and map useful points if needed*/
      if(vox->savePts){
        saveTLSpoints(tempTLS,tInd,vox,scan,xCent,yCent,zCent,map,fInd);
      }

      /*clear voxel intersection arrays*/
      TIDY(rangeList);
      TIDY(voxList);
    }/*beam loop*/

    /*reallocate*/
    if((scan->nPoints>0)&&(scan->nPoints<tempTLS->nPoints)){
      if(!(scan->point=(tlsPoint *)realloc(scan->point,scan->nPoints*sizeof(tlsPoint)))){
        Rprintf("Error in reallocation, allocating %lu\n",scan->nPoints*sizeof(tlsPoint));
        exit(1);
      }
    }else if(scan->nPoints==0)TIDY(scan->point);

    /*determine gap fraction and normalise parameters*/
    for(vPlace=0;vPlace<vox->nVox;vPlace++){

      if(vox->savePts){
        /*gap fraction for scaling point size. Only used when points are saved*/
        for(k=0;k<map->nIn[vPlace];k++){
          pInd=map->mapPoint[vPlace][k];
          if((vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace])>0.0){
            scan->point[pInd].gap=vox->hits[fInd][vPlace]/(vox->hits[fInd][vPlace]+vox->miss[fInd][vPlace]);
          }else scan->point[pInd].gap=1.0;
        }/*point within voxel loop*/
      }/*if we are saving points*/

      /*normalise mean properties*/
      if(vox->hits[fInd][vPlace]>0.0){
        vox->meanRefl[fInd][vPlace]/=vox->hits[fInd][vPlace];
        vox->meanZen[fInd][vPlace]/=vox->hits[fInd][vPlace];
      }else{
        vox->meanRefl[fInd][vPlace]=vox->meanZen[fInd][vPlace]=-1.0;
      }/*mean property normalisation*/
    }/*voxel loop*/
  }/*voxel bound check*/

  /*tidy temporary space*/
  tempTLS=tidyTLScan(tempTLS);
  return(scan);
}/*readOneTLS*/


/*##################################################################*/
/*see if a file is aptx file*/

char checkIfPtx(char *namen)
{
  /*check last three characters*/
  if(!strncasecmp(&namen[strlen(namen)-4],".ptx",4))return(1);
  else                                              return(0);
}/*checkIfPtx*/


/*##################################################################*/
/*tidy multiple TLS structures*/

tlsScan *tidyTLScans(tlsScan *scans,int nScans)
{
  int i=0;

  if(scans){
    for(i=0;i<nScans;i++){
      if(scans[i].beam){
        for(i=0;i<scans[i].nBeams;i++){
          TIDY(scans[i].beam[i].refl);
          TIDY(scans[i].beam[i].r);
        }
        TIDY(scans[i].beam);
      }
      TIDY(scans[i].point);
      TTIDY((void **)scans[i].matrix,4);
      if(scans[i].ipoo){
        fclose(scans[i].ipoo);
        scans[i].ipoo=NULL;
      }
    }
    TIDY(scans);
  }

  return(scans);
}/*tidyTLScans*/


/*##################################################################*/
/*tidy up TLS sctructure*/

tlsScan *tidyTLScan(tlsScan *scan)
{
  uint32_t i=0;

  if(scan){
    if(scan->beam){
      for(i=0;i<scan->nRead;i++){
        TIDY(scan->beam[i].refl);
        TIDY(scan->beam[i].r);
      }
      TIDY(scan->beam);
    }
    TTIDY((void **)scan->matrix,4);
    TIDY(scan->point);
    if(scan->ipoo){
      fclose(scan->ipoo);
      scan->ipoo=NULL;
    }
    TIDY(scan);
  }
  return(scan);
}/*tidyTLSscan*/


/*the end*/
/*##################################################################*/

