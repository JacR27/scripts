#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

#define MAXLINES 10000
#define BLOCK_SIZE 512

void insertionsort(char *v[], int len);

int readlines(char *lineptr[], char *text, char *stor, int32_t l_text);

void writelines(char *lineptr[],int nlines);


main()
{
  //read bam header
  char magic[4];
  int32_t l_text;
  char *text;
  
  fread(magic,sizeof(char),4,stdin);
  fread(&l_text, sizeof(int32_t),1, stdin);
  text = malloc(sizeof(char) * l_text);
  fread(text, sizeof(char),l_text,stdin);
  
  //sort header
  int nlines;
  char *stor;
  char *lineptr[MAXLINES];

  stor = malloc(sizeof(char) * l_text + MAXLINES);
  nlines = readlines(lineptr, text, stor, l_text);
  insertionsort(lineptr,nlines);
  
  //write header
  fwrite(magic,1,4,stdout);
  fwrite(&l_text, sizeof(int32_t),1, stdout);
  writelines(lineptr,nlines);
  free(text);
  free(stor);


  // print the rest of the file to stdout
  char buffer[BLOCK_SIZE];
  for (;;){
    size_t bytes = fread(buffer,  sizeof(char),BLOCK_SIZE,stdin);
    fwrite(buffer, sizeof(char), bytes, stdout);
    fflush(stdout);
    if (bytes < BLOCK_SIZE)
      if (feof(stdin))
	break;
  }
}



void insertionsort(char *v[], int len) //becaues it is stable 
{
  int i, j;
  char *temp;
  for (i = 1; i < len; i++){
    temp = v[i];
    j = i - 1;
    while ((strncmp(temp, v[j], 3) < 0) && (j >= 0)){
      v[j + 1] = v[j];
      j = j - 1;
    }
    v[j + 1] = temp;
  }
}


int readlines(char *lineptr[], char *text, char *stor, int32_t l_text)
{
  int i, nlines;
  int isnewline;
  nlines = 0;
  for (i = 0; i < l_text; i++){
    if (isnewline){
      lineptr[nlines] = &stor[i+nlines];
      isnewline = 0;
    }
    stor[i+nlines] = text[i];
    if (text[i] == '\n'){
      nlines++;
      stor[i+nlines] = '\0';
      isnewline = 1;
    }
  }
  return nlines;
}


void writelines(char *lineptr[],int nlines)
{
  int i;
  for(i=0;i<nlines;i++)
    printf("%s",lineptr[i]);
}



