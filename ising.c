#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <argp.h>
#include "mt19937ar.h"           //mersenne twister random number generator
#include "define.h"

typedef enum boolean {false=0,true=1} boolean;

const char *argp_program_version =
"ising Monte Carlo code v1.0 by Richard Tjornhammar";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[2];            /* ARG1 and ARG2 */
  int verbose;              /* The -v flag */
  char *outfile,*infile;    /* Argument for -o and -i*/
};


typedef struct res{
  float m[(LS+1)*(RD+1)];
} results;

typedef struct lat{
  int* S;
  int* nn;
  int  L;
  float m;
  float e;
} lattice;

typedef struct tem{
  int t;
  float T;
  float P[DIM*2+1];
} temperature;

typedef struct sim{ //for saving data to binary
  results RES;
  lattice LAT;
  temperature TEM;
} simulation;

typedef lattice *platt;
typedef temperature *ptemp;
typedef results *pres;

static struct argp_option options[] =
{
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {"input",  'i', "INFILE", 0,
   "Input from INFILE instead"},
  {"output",  'o', "OUTFILE", 0,
   "Output to OUTFILE instead of to standard output"},
  {0}
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'v':
      arguments->verbose = 1;
      break;
    case 'o':
      arguments->outfile = arg;
      break;
    case 'i':
      arguments->infile = arg;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 2)
	{
	  argp_usage (state);
	}
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static char args_doc[] = "ARG1 ARG2";
static char doc[] =
  "AN ISING MODEL IN C";

static struct argp argp = {options, parse_opt, args_doc, doc};

void ising_error(char* str)
{
  fprintf(stderr,"ERROR:: %s\n",str);
  exit(0);
}

void sort (float *AR, int nrEntries, int whc, int nrc) {
  int sa,sb,s,i,j,k;
  float ft,ft2;

  for(sa = 1; sa < nrEntries; ++sa)
    for(sb = nrEntries-1; sb >= sa; --sb) {
      // compare adjacent elements
      if(AR[(int)(sb - 1)*(int)pow(nrc,whc-1)] > AR[(int)sb*(int)pow(nrc,whc-1) ]) {
	// exchange all elements between rows
	for(i=0;i<nrc;i++){
	  ft = AR[(int)(sb - 1)*(int)pow(nrc,i) ];
	  AR[(int)(sb - 1)*(int)pow(nrc,i) ] = AR[(int)sb*(int)pow(nrc,i) ];
	  AR[(int)sb*(int)pow(nrc,i) ] = ft;
	}
      }
    }
}

void init(platt lat){
  int   i,j,k,n,m,L,N,s3,s2,s1,tmp,D;
  float T;

  L=lat->L;
  D=2*DIM;
  lat->S=malloc(sizeof(int)*pow(L,DIM));
  lat->nn=malloc(sizeof(int)*pow(L,DIM)*D);

  N=(int)pow(L,DIM)-1;
  s1=L;s2=SQ(L);s3=CB(L);

  lat->m=pow(L,DIM);
  lat->e=pow(L,DIM);

  for(i=0;i<L;i++)
    for(j=0;j<L;j++)
      switch(DIM){
      case 3:
	for(k=0;k<L;k++){
	  //fprintf(stdout,"%d %d %d\n",i,j,k);
	  lat->S[i+s1*j+s2*k]=1;//2*((i+s1*j+s2*k)%2)-1;

	  m=(i+s1*j+s2*k)*D;
	  lat->nn[m]  = (i==(L-1))?(0+s1*j+s2*k):(i+s1*j+s2*k+1);
	  lat->nn[m+1]= (i==0)?(L-1+s1*j+s2*k):(i+s1*j+s2*k-1);
	  lat->nn[m+2]= (j==(L-1))?(i+s2*k):(i+s1*j+s2*k+1*s1);
	  lat->nn[m+3]= (j==0)?((L-1)*s1+i+s2*k):(i+s1*j+s2*k-1*s1);
	  lat->nn[m+4]= (k==(L-1))?(i+s1*j):(i+s1*j+s2*k+1*s2);
	  lat->nn[m+5]= (k==0)?((L-1)*s2+i+j*s1):(i+s1*j+s2*k-1*s2);
	}
	break;
      case 2:
	lat->S[i+s1*j]=1;//2*((i+s1*j)%2)-1;

	m=(i+s1*j)*D;
	lat->nn[m]  = (i==(L-1))?(s1*j):(i+s1*j+1);
	lat->nn[m+1]= (i==0)?(L-1+s1*j):(i+s1*j-1);
	lat->nn[m+2]= (j==(L-1))?(i):(i+s1*j+1*s1);
	lat->nn[m+3]= (j==0)?((L-1)*s1+i):(i+s1*j-1*s1);
	break;
      default:
	ising_error("DIMENSION IS WRONG");
      }
}

void temper(ptemp TEMP, int t){           //SETTING PROBABILITIES FOR THIS T
  int i;
  TEMP->t=t;
  TEMP->T=(float)t/(float)NT*2.0+(1.5)*((float)DIM-1.0);
  for(i=0;i<=(DIM*2+1);i++){
    TEMP->P[i]=exp(-2*i/(TEMP->T));
  }
}

void sweep(platt lat, ptemp temp){
  int i,j,k,l;
  int s,w;
  int N,L,NSW,D;

  L=lat->L;
  l=L-1;
  D=2*DIM;

  switch(DIM){
  case 3:
    N=CB(L)-1;
    for(i=0;i<N;i++){
      j=(int)(rnd()*N);//(lat->L)*t

      s=lat->S[j];
      w=s*( lat->S[lat->nn[j*D]]+lat->S[lat->nn[j*D+1]]+lat->S[lat->nn[j*D+2]]+lat->S[lat->nn[j*D+3]]+lat->S[lat->nn[j*D+4]]+lat->S[lat->nn[j*D+5]]);
      if(w<=0 || temp->P[w]>rnd()) {
	lat->S[j]=-s;
	lat->m+=2*lat->S[j];
	lat->e+=2*w;
      }
    }
    break;
  case 2:
    N=SQ(L)-1;
    for(i=0;i<N;i++){
      j=(int)(rnd()*N);

      s=lat->S[j];
      w=s*( lat->S[lat->nn[j*D]]+lat->S[lat->nn[j*D+1]]+lat->S[lat->nn[j*D+2]]+lat->S[lat->nn[j*D+3]]);
      if(w<=0 || temp->P[w]>rnd()) {
	lat->S[j]=-s;
	lat->m+=2*lat->S[j];
	lat->e+=2*w;
      }
    }
    break;
  default:
    break;
  }
}

void equi(platt lat, ptemp temp){        //equilibration of the lattice
  int   i,j,k,L,N;
  float T;

  for(i=0;i<NWARM;i++)
    sweep(lat,temp);

}

void prod(platt lat, ptemp temp, pres res){  //excercise the lattice
  int   i,j,k,L,N;
  float m=0,ma=0,ma2=0,ma4=0,e=0,e2=0;
  L=lat->L;
  N=(DIM==3)?CB(L):SQ(L);

  for(i=0;i<NRUN;i++){
    sweep(lat,temp);
    m+=lat->m;
    ma+=abs(lat->m);
    ma2+=SQ(lat->m);
    ma4+=SQ(SQ(lat->m));
    e+=lat->e;
    e2+=SQ(lat->e);
  }
  m/=NRUN;
  ma/=NRUN;
  ma2/=NRUN;
  ma4/=NRUN;
  e/=NRUN;
  e2/=NRUN;

  j=temp->t;
  //fprintf(stdout,"\n%d\n\n",j);
  res->m[j*RD]=temp->T;
  res->m[j*RD+1]=(float)m/(float)N;
  res->m[j*RD+2]=(float)ma/(float)N;
  res->m[j*RD+3]=(float)(ma2-SQ(ma))/(float)N/temp->T;     //chi
  res->m[j*RD+4]=(float)e/(float)N;                      //e
  res->m[j*RD+5]=1.0-pow( (float)ma4/(float)ma2,0.66666);//binder
}

void coll(platt lat,pres res, FILE *stream){         //sort data and print
  int i,j,k;

  //sort
  for(i=0;i<NT;i++){
    fprintf(stream,"%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",res->m[i*RD],res->m[i*RD+1],res->m[i*RD+2],res->m[i*RD+3],res->m[i*RD+4],res->m[i*RD+5]);
  }

}

void printLattice(platt lat, FILE *stream){
  //plots the first slice of any lattice
  int i,j,l;

  l=lat->L;
  fprintf(stream,"\n\n(%d,%d)-LATTICE >> \n",l,l);
  for(i=0;i<lat->L;i++){
    for(j=0;j<lat->L;j++)
      fprintf(stream,"% 2d ",lat->S[i+j*l]);
    fprintf(stream,"\n");
  }
}

void printNN(platt lat, FILE *stream)
{
  int N,i,t,D;
  D=2*DIM;

  N=pow(lat->L,DIM);
  for(t=0;t<N;t++){
    fprintf(stream,"\n%d:: ",t);
    for(i=0;i<D;i++)
      fprintf(stream," %d ",lat->nn[D*t+i]);
  }
  fprintf(stream,"\n\n");
}

void kill(platt lat, pres res, ptemp temp){
  //free(lat->S);
  //free(lat->nn);
  //  free(res);
  //  free(temp);
}

int main(int argc, char **argv)
{
  //IO
  struct arguments arguments;
  FILE *outstream;
  int threads;

  //CONTAINERS
  platt lat;
  pres res;
  ptemp tem;

  //RUN THINGS
  int   i,t,N;

  lat = malloc(sizeof(lattice));
  res = malloc(sizeof(results));
  tem = malloc(sizeof(temperature));

  arguments.outfile = NULL;
  arguments.verbose = 0;

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  if (arguments.outfile)
    outstream = fopen (arguments.outfile, "w");
  else
    outstream = stdout;

  lat->L = atoi(arguments.args[0]);
  if(!(lat->L<LS))
    ising_error("LATTICE IS TOO LARGE");

  threads=atoi(arguments.args[1]);

  init(lat);
  if(arguments.verbose){
    printLattice(lat,outstream);
    printNN(lat,outstream);
  }

  for(t=0;t<NT;t++){         //Temperature loop
    temper(tem,t);
    equi(lat,tem);           //equilibration of the lattice
    prod(lat,tem,res);       //excercise the lattice
  }
  coll(lat,res,outstream);   //sort data and print
  if(arguments.verbose)
    printLattice(lat,outstream);

  //if (arguments.outfile)
    //  fclose(outstream);
  //kill(lat,res,tem);

  return 0;
}
