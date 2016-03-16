#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <string.h>


#define MAXNAMELEN 10000
#define MAXLINES 10000
#define MAXSTRINGLEN 100
#define MAXTAGS 100
#define MAXLEN 1000
#define BLOCK_SIZE 512
struct tag readtag();

void printtag(struct tag *t1);

void getstring(char *);

int compare(const void *s1, const void *s2);

void printtag(struct tag *t1);

void insertion(char *v[], int len);

int readlines(char *lineptr[], char *text, char *stor, int32_t l_text);

void writelines(char *lineptr[],int nlines);

void outputtag(struct tag *t1);

union Uvalue {
  int32_t ival;
  char sval[MAXSTRINGLEN];
};

struct tag {
  int len;
  int new_len;
  int order;
  char tag_name[3];
  char tag_value_type;
  char old_tag_value_type;
  union Uvalue value;
};

main()
{
  //read bam header
  int ref;
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
  insertion(lineptr,nlines);
  
  
  //write header
  fwrite(magic,1,4,stdout);
  fwrite(&l_text, sizeof(int32_t),1, stdout);
  writelines(lineptr,nlines);
  free(text);
  free(stor);

  char buffer[BLOCK_SIZE];
  while (1){
    size_t bytes = fread(buffer,  sizeof(char),BLOCK_SIZE,stdin);
    fwrite(buffer, sizeof(char), bytes, stdout);
    fflush(stdout);
    if (bytes < BLOCK_SIZE)
      if (feof(stdin))
	break;
  }
}

parsefile()
{

  //read bam header
  int ref;
  char magic[4];
  int32_t l_text;
  char *text;
  int32_t n_ref;
  int32_t l_name;
  char *name;
  int32_t l_ref;
  
  fread(magic,sizeof(char),4,stdin);
  fread(&l_text, sizeof(int32_t),1, stdin);
  text = malloc(sizeof(char) * l_text);
  fread(text, sizeof(char),l_text,stdin);
  fread(&n_ref, sizeof(int32_t),1, stdin);
  
  
  //sort header
  int nlines;
  char *stor;
  char *lineptr[MAXLINES];

  stor = malloc(sizeof(char) * l_text + MAXLINES);
  nlines = readlines(lineptr, text, stor, l_text);
  insertion(lineptr,nlines);
  
  
  //write header
  fwrite(magic,1,4,stdout);
  fwrite(&l_text, sizeof(int32_t),1, stdout);
  writelines(lineptr,nlines);
  fwrite(&n_ref, sizeof(int32_t),1, stdout);
  free(text);
  free(stor);
  
  //read and write references
  char *refs;
  refs = malloc(n_ref * (sizeof(int32_t)*2 + MAXNAMELEN));
  for (ref=0; ref<n_ref; ++ref){
    fread(&l_name, sizeof(int32_t),1, stdin);
    name = malloc(sizeof(char) * l_name);
    fread(name, sizeof(char), l_name,stdin);
    fread(&l_ref, sizeof(int32_t), 1, stdin);
    
    fwrite(&l_name, sizeof(int32_t),1, stdout);
    fwrite(name, sizeof(char), l_name,stdout);
    fwrite(&l_ref, sizeof(int32_t), 1, stdout);
    free(name);
  }
  free(refs);
  

  //read bam record 
  int32_t block_size;
  int32_t refID;
  int32_t pos;
  uint32_t bin_mq_nl;
  uint32_t flag_nc;
  int32_t l_seq;
  int32_t next_refID;
  int32_t next_pos;
  int32_t tlen;
  char *read_name;
  uint32_t *cigar;
  uint8_t *seq;
  char *qual;

  char l_read_name;
  int n_cigar_op;
  int block_size_no_tags;
  int bytes_left_in_block;
  int seq_size;

  int32_t new_block_size;
  int count;
  struct tag new1_tags[MAXTAGS];
  int new_count;
  int n_tags;
  int i;
  
  while (1){
    fread(&block_size, sizeof(int32_t), 1, stdin);
    if(feof(stdin))
      break;
    fread(&refID, sizeof(int32_t), 1, stdin);
    fread(&pos, sizeof(int32_t), 1, stdin);
    fread(&bin_mq_nl, sizeof(uint32_t), 1, stdin);
    fread(&flag_nc, sizeof(uint32_t), 1, stdin);
    fread(&l_seq, sizeof(int32_t), 1, stdin);
    fread(&next_refID, sizeof(int32_t), 1, stdin);
    fread(&next_pos, sizeof(int32_t), 1, stdin);
    fread(&tlen, sizeof(int32_t), 1, stdin);
    
    l_read_name = bin_mq_nl%256;
    n_cigar_op = flag_nc%65536;
    seq_size = (l_seq+1)/2;
    
    read_name = malloc(sizeof(char) * l_read_name);
    fread(read_name, sizeof(char), l_read_name, stdin);
    
    cigar = malloc(sizeof(uint32_t) * n_cigar_op);
    fread(cigar, sizeof(uint32_t), n_cigar_op, stdin);
  
    seq = malloc(sizeof(uint8_t)*seq_size);
    fread(seq, sizeof(uint8_t), seq_size, stdin);
    
    qual = malloc(sizeof(char) * l_seq);
    fread(qual, sizeof(char), l_seq, stdin);
    
    //read and sort tags and modify tags
    block_size_no_tags = 32 + l_read_name + 4*n_cigar_op + seq_size + l_seq;
    bytes_left_in_block = block_size - block_size_no_tags;
    count = 0;
    new_count = 0;
    n_tags = 0;
  
    while (bytes_left_in_block != count){
      new1_tags[n_tags] = readtag();
      count += new1_tags[n_tags].len;
      new_count += new1_tags[n_tags].new_len;
      ++n_tags;
  }
    
    new_block_size = block_size_no_tags + new_count; 
    qsort(new1_tags, n_tags, sizeof(struct tag), compare);
   
    //write out record and tags
    fwrite(&new_block_size, sizeof(int32_t), 1, stdout);
    fwrite(&refID, sizeof(int32_t), 1, stdout);
    fwrite(&pos, sizeof(int32_t), 1, stdout);
    fwrite(&bin_mq_nl, sizeof(uint32_t), 1, stdout);
    fwrite(&flag_nc, sizeof(uint32_t), 1, stdout);
    fwrite(&l_seq, sizeof(int32_t), 1, stdout);
    fwrite(&next_refID, sizeof(int32_t), 1, stdout);
    fwrite(&next_pos, sizeof(int32_t), 1, stdout);
    fwrite(&tlen, sizeof(int32_t), 1, stdout);
    fwrite(read_name, sizeof(char), l_read_name, stdout);
    fwrite(cigar, sizeof(uint32_t), n_cigar_op, stdout);
    fwrite(seq, sizeof(uint8_t), seq_size, stdout);
    fwrite(qual, sizeof(char), l_seq, stdout);    
  
    for (i = 0; i<n_tags; i++){
      outputtag(new1_tags+i);
    }
  
    free(read_name);
    free(cigar);
    free(seq);
    free(qual);
  }
}


struct tag readtag()
{
  struct tag tg;
  int len = 3;
  tg.new_len = 3;
  char tag_name[3] = "\0\0";
  char tag_value_type;
  union Uvalue value;
  int tmp= 0;
  fread(tag_name, sizeof(char), 2, stdin);
  if (strcmp(tag_name,"SM") == 0)
    tg.order =1;
  if (strcmp(tag_name,"AS") == 0)
    tg.order =2;
  if (strcmp(tag_name,"RG") == 0)
    tg.order =3;
  if (strcmp(tag_name,"NM") == 0)
    tg.order =4;
  if (strcmp(tag_name,"BC") == 0)
    tg.order =5;
  if (strcmp(tag_name,"OC") == 0)
    tg.order =6;
  if (strcmp(tag_name,"SA") == 0)
    tg.order =7;
  fread(&tag_value_type, sizeof(char), 1, stdin);
  switch (tag_value_type){
  case 'Z':
    getstring(value.sval);
    len += (strlen(value.sval)+1);
    tg.new_len += (strlen(value.sval)+1);
    tg.old_tag_value_type = tag_value_type;
    break;
  case 'c':
    fread(&tmp, sizeof(int8_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int8_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'C':
    fread(&tmp, sizeof(uint8_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint8_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 's':
    fread(&tmp, sizeof(int16_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int16_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'S':
    fread(&tmp, sizeof(uint16_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint16_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'i':
    fread(&tmp, sizeof(int32_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(int32_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  case 'I':
    fread(&tmp, sizeof(uint32_t), 1, stdin);
    value.ival = tmp;
    len += sizeof(uint32_t);
    tg.old_tag_value_type = tag_value_type;
    tag_value_type = 'i';
    tg.new_len += sizeof(int32_t);
    break;
  }
  strcpy(tg.tag_name,tag_name);
  tg.tag_value_type = tag_value_type;
  tg.value = value;
  tg.len = len;
  return tg;
}  

void printtag(struct tag *t1)
{
  fprintf(stderr, "tag len: %d\n ",t1->len);
  fprintf(stderr, "new tag len: %d\n ",t1->new_len);
  fprintf(stderr, "tag name: ");
  fwrite(t1->tag_name, sizeof(char),2,stderr);
  fprintf(stderr, "\n");
  fprintf(stderr, "tag order %d\n", t1->order);
  fprintf(stderr, "tag value type: %c\n ",t1->tag_value_type);
  fprintf(stderr, "old tag value type: %c\n ",t1->old_tag_value_type);
  
  if ((t1->tag_value_type) == 'Z')
    fprintf(stderr, "tag value: %s\n ",(t1->value).sval);
  else
    fprintf(stderr, "tag value: %d\n ",(t1->value).ival);
  fprintf(stderr, "\n");
}

void outputtag(struct tag *t1)
{
  fwrite(t1->tag_name, sizeof(char),2,stdout);
  fprintf(stdout, "%c",t1->tag_value_type);
  if ((t1->tag_value_type) == 'Z')
    fwrite((t1->value).sval,sizeof(char),(t1->len)-3,stdout);
  else{
    int32_t i = (t1->value).ival + 0;
    fwrite(&i,sizeof(int32_t),1,stdout);
  }
}

void getstring(char *s)
{
  char c;
  while (1){
    c = getchar();
    *s++ = c;
    if (c == '\0')
      break;
  }
}


int compare(const void *s1, const void *s2)
{
  struct tag *t1 = (struct tag *)s1;
  struct tag *t2 = (struct tag *)s2;
  int tagcompare = t1->order - t2->order;
  return tagcompare;
}


void insertion(char *v[], int len)
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
    //if (strstr(lineptr[i], "scramble") != NULL)
    printf("%s",lineptr[i]);
}



