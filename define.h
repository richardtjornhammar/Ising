//compile with gcc *.c -lm -o ising
//             gcc mt19937ar.* ising.c -lm -o ising
//             gcc -fopenmp mt19937ar.* ising_omp.c -lm -o ising
#define rnd	genrand_real1
#define PI	3.1415926535898
#define DIM	3
#define RD      6
#define NWARM   1000
#define NRUN    1000
#define	SQ(X)	(X)*(X)
#define CB(X)   (X)*(X)*(X)
#define SGN(X)  X>0?1:-1
#define MAX	1000
#define LS      100
#define NT      100
