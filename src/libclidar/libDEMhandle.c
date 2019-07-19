#include <Rinternals.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "stdint.h"
#include "tools.h"
#include "libDEMhandle.h"


/*#########################*/
/*# Functions to handle  #*/
/*# DEMs to process lidar #*/
/*#########################*/


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



/*#######################################*/
/*read an ASCII DEM*/

demStruct *readAscDEM(char *namen,double minX,double minY,double maxX,double maxY)
{
  int i=0,j=0,place=0;
  int maxLen=0,lineN=0;
  demStruct *dem=NULL;
  char *line=NULL,temp1[100],temp2[100];
  char *token=NULL;
  FILE *ipoo=NULL;

  maxLen=100000;
  line=challoc((uint64_t)maxLen,"line",0);

  if(!(dem=(demStruct *)calloc(1,sizeof(demStruct)))){
    Rprintf("error demStruct allocation.\n");
    exit(1);
  }
  if((ipoo=fopen(namen,"r"))==NULL){
    Rprintf("Error opening dem file %s\n",namen);
    exit(1);
  }

  lineN=0;
  while(fgets(line,maxLen,ipoo)!=NULL){
    if(lineN<6){ /*read the header*/
      if(sscanf(line,"%s %s",temp1,temp2)==2){
        if(!strncasecmp(line,"ncols",5))dem->nX=atoi(temp2);
        else if(!strncasecmp(line,"nrows",5))dem->nY=atoi(temp2);
        else if(!strncasecmp(line,"xllcorner",9))dem->minX=atof(temp2);
        else if(!strncasecmp(line,"yllcorner",9))dem->minY=atof(temp2);
        else if(!strncasecmp(line,"cellsize",8))dem->res=atof(temp2);
        else if(!strncasecmp(line,"NODATA_value",12)){
          dem->noData=atof(temp2);
          j=dem->nY-1;
          dem->z=dalloc(dem->nX*dem->nY,"dem",0);
        }
      }
    }else{  /*read data*/
      token=strtok(line," ");
      i=0;
      while(token!=NULL) {
        place=j*dem->nY+i;
        dem->z[place]atof(token);
        i++;
      }
      j--;
    }
    lineN++;
  }

  /*allocate space*/


  TIDY(line);
  if(ipoo){
    fclose(ipoo);
    ipoo=NULL;
  }
  return(dem);
}/*readAscDEM*/

/*the end*/
/*#######################################*/

