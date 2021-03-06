#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <zlib.h>
#include <ctype.h>
#include <unistd.h>
#include "deflate.h"

int mode, lvalue, tvalue, wvalue, svalue, cvalue, nFiles, cfvalue, level, strat;
int fflag, rflag, mflag, jflag, iflag, dflag;
char bclFormat[100];
char outExtention[100];
char printFile[50];
main( int argc, char *argv[])

    {
      extern int mode, lvalue, tvalue, wvalue, svalue, cvalue, level, strat;
      extern char bclFormat[];
      

      level = 6;
      strat = Z_HUFFMAN_ONLY;
      lvalue = 1;
      tvalue = 24;
      wvalue = 2;
      svalue = 2;
      cvalue = 302;
      cfvalue = 1;
      mode = 0;
      strcpy(bclFormat, "bcl");
      strcpy(outExtention, "bcl.6hf.gz");

      extern int fflag, rflag, mflag, jflag, iflag, dflag;
      fflag = 0; //filtering
      rflag = 0; // reduceingResolution
      mflag = 0; // remapping
      jflag = 0; // joining
      iflag = 0; // spliting
      dflag = 0; // demultiplexing

      
      int index;
      int c;

      opterr = 0;

      while ((c = getopt (argc, argv, "p:o:l:t:w:s:c:z:frmjidg:a:h:")) != -1)
	switch (c)
	  {
	  case 'o':
	    mode = atoi(optarg);
	    if ((mode == 1) || (mode == 3))
	      fflag = 1;
	    break;
	  case 'l':
	    lvalue = atoi(optarg);
	    break;
	  case 't':
	    tvalue = atoi(optarg);
	    break;
	  case 'w':
	    wvalue = atoi(optarg);
	    break;
	  case 's':
	    svalue = atoi(optarg);
	    break;
	  case 'c':
	    cvalue = atoi(optarg);
	    break;
	  case 'a':
            level = atoi(optarg);
            break;
	  case 'g':
            strat = atoi(optarg);
            break;
	  case 'p':
	    strcpy(printFile, optarg);
	  case 'z':
	    strcpy(bclFormat, optarg);
	    printf("%s\n",bclFormat);
	    break;
	  case 'h':
	    strcpy(outExtention, optarg);
	    printf("%s\n",bclFormat);
	  break;
	  case 'f':
	    fflag = 1;
	    break;
	  case 'r':
	    rflag = 1;
	    break;
	  case 'm':
	    mflag = 1;
	    break;
	  case 'j':
	    jflag = 1;
	    break;
	  case 'i':
	    iflag = 1;
	    break;
	  case 'd':
	    dflag = 1;
	    break;
	  case '?':
	    if (optopt == 'l' || optopt == 'o' || optopt == 't' || optopt == 's' || optopt == 'f' || optopt == 'c')
	      fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	    else if (isprint (optopt))
	      fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	    else
	      fprintf (stderr,
		       "Unknown option character `\\x%x'.\n",
		       optopt);
	    return 1;
	  default:
	    abort ();
	  }

      printf ("lane = %d, cycles = %d, ntiles = %d, nswaths = %d, nsurfaces = %d, bclformat = %s, filter = %d, reduce resolution = %d, remap = %d, join = %d, split = %d, demultiplexed = %d\n",
	      lvalue ,cvalue, tvalue , wvalue , svalue ,bclFormat, fflag , rflag , mflag, jflag, iflag, dflag);

      for (index = optind; index < argc; index++)
	printf ("Non-option argument %s\n", argv[index]);
      
      	
      frr();
      
      return 0;
    }



void frr(void)
{
  #define MAXFOLDERNAMELEN 8
  printf("starting main\n");
  extern int cvalue;
  //extern char **folderNames;
  int nCycles;
  nCycles = (cvalue + 1) - cfvalue;
  printf("%d\n", nCycles);
  char folderNames[nCycles][MAXFOLDERNAMELEN];
  int f;
  int t;
  for (f = cfvalue, t = 0; f <= cvalue; ++f, ++t){
    printf("%d\n",f);
    sprintf(folderNames[t],"C%d.1/",f);
  }



#define MAXFILENAMELEN 50


  int i, j, k , n;
  nFiles= (svalue * wvalue * tvalue);
  char bclNames[nFiles][MAXFILENAMELEN];
  extern int svalue, tvalue, wvalue, lvalue, nFiles;  
  n = 0;
  printf("nfiles %d\n",nFiles);
  for (i = 1; i <= svalue ; ++i){
    for (j = 1; j <= wvalue; ++j){
      for (k = 1; k <= tvalue; ++k){
	snprintf(bclNames[n],MAXFILENAMELEN,"s_%d_%d%d%02d.",lvalue,i,j,k);
	n++;
      }
    }
  }
  
  printf("got file names");

  

  #define MAXEXTENTIONLEN 100
    
  
  int u;
  int w;
  char inpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+MAXEXTENTIONLEN];
  char outpath[MAXFILENAMELEN+MAXFOLDERNAMELEN+MAXEXTENTIONLEN];
  for (u=0;u<nCycles;u++){
    for( w=0;w<nFiles;w++){
      strcpy(inpath,folderNames[u]);
      strcat(inpath,bclNames[w]);
      strcat(inpath,bclFormat);
      printf("%s\n", inpath);
      strcpy(outpath,folderNames[u]);
      strcat(outpath,bclNames[w]);
      strcat(outpath,outExtention);
      printf("%s\n", outpath);
      
      FILE *inputFile;
      FILE *outputFile;
      outputFile = fopen(outpath,"w");
      inputFile = fopen(inpath,"r");
      def(inputFile,outputFile,level,strat);

      fclose(inputFile);
      fclose(outputFile);




    }
  }
}



