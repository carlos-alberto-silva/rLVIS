int compFloat(const void *x,const void *y)
{
  int returning=0;

  if((*(float *)x)<(*(float *)y))     returning=-1;
  else if((*(float *)x)>(*(float *)y))returning=1;
  else                                returning=0;

  return(returning);
}/*compFloat*/



/*#########################################################################*/
/*allocate a double arrays, checking for errors*/

double *dalloc(int length,char *namen,int n)
{
  double *jimlad=NULL;
  if(!(jimlad=(double *)calloc(length,sizeof(double)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %d\n",length);
    exit(1);
  }
  return(jimlad);
}/*dalloc*/



float *falloc(uint64_t length,char *namen,int n)
{
  float *jimlad=NULL;
  if(!(jimlad=(float *)calloc(length,sizeof(float)))){
    fprintf(stderr,"error in %s array allocation %d\n",namen,n);
    fprintf(stderr,"allocating %lu\n",length);
    exit(1);
  }
  return(jimlad);
}/*falloc*/


float *medianFloat(float *jimlad,int width,int length)
{
  int i=0,j=0,place=0;
  int halfLength=0,nIn=0;
  float *filtered=NULL,*temp=NULL;
  int compFloat(const void *x,const void *y);  /*function needed by qsort()*/

  if(width<3){
    fprintf(stderr,"Median filter won't work with a width of %d for float\n",width);
    exit(1);
  }

  filtered=falloc((uint64_t)length,"median filter",0);
  temp=falloc((uint64_t)width,"temporary median",0);

  halfLength=width/2;
  for(i=0;i<length;i++){
    nIn=0;
    for(j=0;j<width;j++){
      place=i+j-halfLength;
      if((place>0)&&(place<length)){
        temp[nIn++]=jimlad[place];
      }
    }/*median loop*/
    for(j=nIn;j<width;j++)temp[j]=9999.0;    /*pad end if needed so that they get put to the back during re-ordering*/
    qsort(temp,width,sizeof(float),compFloat);  /*put the contents of temp in order*/
    if(nIn>0){
      filtered[i]=temp[(int)(nIn/2)];
    }else{
      fprintf(stderr,"That didn't work, number in %d for bin %d of %d\n",nIn,i,length);
      exit(1);
    }
  }  /*array loop*/

  TIDY(temp);
  return(filtered);
}/*medianFloat*/
