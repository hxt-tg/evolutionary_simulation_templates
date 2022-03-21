#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <error.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#if defined(__linux__)
#include <sys/stat.h>
#include <linux/limits.h>
#elif defined(__MINGW32__) || defined(__MINGW64__)
#include <sys/stat.h>
#include <limits.h>
#else              // defined(_WIN32) || defined(_WIN64)
#include <limits.h>
#endif

/* ====================== Settings ======================= */
/* Global settings (Comment to disable) */
//#define USE_MATRIX_INIT             /* Call init_by_st_mat() function to initialize strategies (configure: ST_MAT) */
//#define RUN_SPEC_PARAM              /* Running under special parameters, uncomment for loop */
#define REC_STEPS                   /* Record data per step */
#define REC_SNAPSHOT                /* Record snapshots */
/* Global settings */
#define N_PARAM       1             /* Number of Looping variables */
#define RANDSEED      3145215       /* Fixed random seed, use time if -1 */
/* Network settings */
#define L             100           /* Side length of a network */
#define SIZE          (L * L)       /* Total number of nodes */
#define N_NEIGH       4             /* Number of neighbors for every node */
/* Monte Carlo settings */
#define K             0.1           /* Fermi coefficient */
#define MC_STEPS      10000         /* Number of MC steps */
#define AVG_STEPS     5000          /* Number of steps at the end for averaging */
#define REFRESH_FREQ  100           /* Control the frequency of refresh screen, 1 for all the time */
#define TRY_TIMES     1             /* Average result in tried times */
/* Running parameters */
double b_start = 1, b_end = 1.5, b_interval = 0.05;
/* =================== End of settings =================== */

/* Running under special parameters. (Enabled by RUN_SPEC_PARAM) */
double SPEC_PARAM_LIST[][N_PARAM] = {
    /* {b}, */
    {1.00},
    {1.02},
    {1.06},
    {1.20},
};


/* =================== Random function =================== */
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df          /* constant vector a */
#define UPPER_MASK 0x80000000        /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff        /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define NAMEOUT     "K4b075r5Q2"
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)
static unsigned long mt[NN];        /* the array for the state vector  */
static int mti = NN + 1;            /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed) {
    int i;
    for (i = 0; i < NN; i++) {
        mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
        mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
    }
    mti = NN;
}
void lsgenrand(unsigned long seed_array[]) {
    int i; for (i = 0; i < NN; i++) mt[i] = seed_array[i]; mti = NN;
}
double genrand() {
    unsigned long y;
    static unsigned long mag01[2] = { 0x0, MATRIX_A };
    if (mti >= NN) {
        int kk;
        if (mti == NN + 1) sgenrand(4357);
        for (kk = 0; kk < NN - MM; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (; kk < NN - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MM - NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[NN - 1] = mt[MM - 1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;
}
double randf() { return ((double)genrand() * 2.3283064370807974e-10); }
long randi(unsigned long LIM) { return((unsigned long)genrand() % LIM); }
/* ==================== End of random ==================== */


/* ================== Utility functions ================== */
/* Buffers and global vars */
char _file_path[PATH_MAX];
unsigned long _rand;

#define _fprint_str(fp, x, interpret) \
    fprintf(fp, "%s,\"%s\",\"(%s)\"\n", #x, x, interpret)
#define _fprint_int(fp, x, interpret) \
    fprintf(fp, "%s,%d,\"(%s)\"\n", #x, x, interpret)
#define _fprint_dbl(fp, x, interpret) \
    fprintf(fp, "%s,%g,\"(%s)\"\n", #x, x, interpret)

void save_meta() {
    FILE *fp = fopen("output/meta.csv", "w");
    _fprint_int(fp, L, "Side length of, square lattice");
    _fprint_int(fp, SIZE, "Number of agents on network");
    _fprint_int(fp, N_NEIGH, "Number of neighbors for every agent");
    _fprint_dbl(fp, K, "Fermi coefficient");
    _fprint_int(fp, (int)_rand, "Random seed");
    _fprint_int(fp, MC_STEPS, "Monte Carlo steps");
    _fprint_int(fp, AVG_STEPS, "Recorded MC steps at the end");
    _fprint_int(fp, TRY_TIMES, "Try times per parameter");
#ifndef RUN_SPEC_PARAM
    fprintf(fp, "b=%g~%g with interval %g.\n", b_start, b_end, b_interval);
#endif
    fclose(fp);
}

void make_dir(const char *path) {
    int status;

#if defined(__linux__)
    status = mkdir(path, 0776);
#elif defined(__MINGW32__) || defined(__MINGW64__)
    status = mkdir(path);
#else
    fprintf(stderr, "Unsupport platform. (Used in make_dir)\n");
    exit(ENOTSUP);
#endif

    if (status && errno != EEXIST) {
        fprintf(stderr, "Cannot mkdir 'output'.\n");
        exit(errno);
    }
}
/* ============== End of utility functions =============== */








double b;
int net[SIZE][N_NEIGH];   /* Players neighbor relations */

/* Properties counter */
typedef enum {
    _C,
    _D,
    N_STRA
} Strategy;

Strategy stra[SIZE];
int stra_cnt[N_STRA] = {0};

/* Payoff matrix and its update (N_STRA*N_STRA matrix) */
double payoff_matrix[N_STRA][N_STRA] = {
    {1, 0},
    {0, 0}
};
/* Call update_matrix(b) after loop for b */
#define update_matrix(b) do {    \
    payoff_matrix[1][0]=b;       \
} while (0)

void prod_neighbors() {
    int i, j, x;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++){
            x = i*L + j;
            net[x][0] = i*L + ((j-1+L)%L);  /* left */
            net[x][1] = ((i-1+L)%L)*L + j;  /* up */
            net[x][2] = ((i+1)%L)*L + j;    /* down */
            net[x][3] = i*L + ((j+1)%L);    /* right */
        }
}

void init_randomly() {
    int i;
    for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
    for (i = 0; i < SIZE; i++) stra_cnt[stra[i] = (Strategy)randi(N_STRA)]++;
}

/* ST_MAT: Configure strategies matrix for initialization. */
#define ST_MAT_ROW  2
#define ST_MAT_COL  2
const char ST_MAT[ST_MAT_ROW][ST_MAT_COL] = {
    { _C, _D},
    { _D, _C},
};

void init_by_st_mat() {
    int i, j, block_row = L/ST_MAT_ROW, block_col = L/ST_MAT_COL;
    for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++) 
            stra_cnt[stra[i*L+j] = ST_MAT[i/block_row][j/block_col]]++;
}

double payoff(int x) {
    double pay = 0;
    int i;
    for (i = 0; i < N_NEIGH; i++)
        pay += payoff_matrix[stra[x]][stra[net[x][i]]];
    return pay;
}

void update_stra(int x) {
    int y = net[x][(int)randi(N_NEIGH)];
    /* Update strategy, if same then exit */
    if (stra[x] == stra[y]) return ;
    if (randf() < 1/(1+exp((payoff(x) - payoff(y))/K))){
        stra_cnt[stra[x]]--;
        stra[x] = stra[y];
        stra_cnt[stra[x]]++;
    }
}

/* Snapshot function temporarily not used */
void snapshot(int try_time, int step) {
#ifndef REC_SNAPSHOT
    return ;
#endif
    if (TRY_TIMES > 1) sprintf(_file_path, "output/snapshots/b=%.2lf_step=%d[%d].csv", b, step+1, try_time);
    else sprintf(_file_path, "output/snapshots/b=%.2lf_step=%d.csv", b, step+1);
    FILE *fp = fopen(_file_path, "w");
    int x;
    for (x = 0; x < SIZE; x++)
        fprintf(fp, "%d%c", stra[x], (x+1)%L?',':'\n');
    fclose(fp);
}

void run(FILE *fp_aver) {
    FILE *fp_steps = NULL;
#ifdef REC_STEPS
    sprintf(_file_path, "output/steps/b=%.2lf.csv", b);
    fp_steps = fopen(_file_path, "w");
    if (TRY_TIMES > 1) fprintf(fp_steps, "try,");
    fprintf(fp_steps, "step,f_C,f_D\n");
#endif
    int i, x, step, rec, ta[N_STRA]={0}, try_time;      /* total amount for stras in all REC steps */
    double f[N_STRA]={0}, af[TRY_TIMES][N_STRA]={{0}};    /* freq and avg_freq */
    update_matrix(b);
    for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
        for (i = 0; i < N_STRA; i++) ta[i] = 0;
        rec = 0;
#ifdef USE_MATRIX_INIT
        init_by_st_mat();
#else
        init_randomly();
#endif
        snapshot(try_time, -1);
        for (step = 0; step < MC_STEPS; step++) {
            for (x = 0; x < SIZE; x++)
                update_stra((int)randi(SIZE));

            for (i = 0; i < N_STRA; i++)
                f[i] = stra_cnt[i]/(double)SIZE;

            if (step > MC_STEPS-AVG_STEPS-1) {
                for (i = 0; i < N_STRA; i++)
                    ta[i] += stra_cnt[i];
                rec++;
            }

            if (TRY_TIMES > 1) fprintf(fp_steps, "%d,", try_time);
            fprintf(fp_steps, "%d,%lf,%lf\n", step, f[_C], f[_D]);

            if (REFRESH_FREQ && step % REFRESH_FREQ == 0) {
                if (TRY_TIMES > 1) printf("\r [%d]", try_time);
                else printf("\r");
                printf(" Step: %6d  f: [%6.2lf%%,%6.2lf%%]", step, f[_C]*100, f[_D]*100);
            }

            if (stra_cnt[_C] == SIZE || stra_cnt[_D] == SIZE) {
                if (step++ < MC_STEPS-AVG_STEPS)
                    for (i = 0; i < N_STRA; i++)
                        af[try_time][i] = stra_cnt[i] ? 1 : 0;
                break;
            }
        }
        if (rec)
            for (i = 0; i < N_STRA; i++)
                af[try_time][i] = (double)ta[i]/(SIZE * rec);

        if (TRY_TIMES > 1) printf("\r [%d]", try_time);
        else printf("\r");
        printf(" b: %.2lf   avg: [%6.2lf%%,%6.2lf%%]       \n", b,
                af[try_time][_C]*100, af[try_time][_D]*100);

    }
#ifdef REC_STEPS
    fclose(fp_steps);
#endif

    double aaf[N_STRA] = {0};       /* Average in all tries */
    for (i = 0; i < N_STRA; i++) {
        for (x = 0; x < TRY_TIMES; x ++)
            aaf[i] += af[x][i];
        aaf[i] /= TRY_TIMES;
    }
    fprintf(fp_aver, "%.2lf,%lf,%lf\n", b, aaf[_C], aaf[_D]);
    if (TRY_TIMES > 1) printf("[%d tries] b: %g  avg: [%6.2lf%%,%6.2lf%%]\n\n",
            TRY_TIMES, b, aaf[_C]*100, aaf[_D]*100);
}

int main() {
    make_dir("output");
    
#ifdef REC_STEPS
    make_dir("output/steps");
#endif

#ifdef REC_SNAPSHOT
    make_dir("output/snapshots");
#endif

    /* Flush file */
    FILE *fp_aver = fopen("output/average.csv", "w");
    fprintf(fp_aver, "b,f_C,f_D\n");

    _rand = RANDSEED == -1 ? time(NULL) : RANDSEED;
    sgenrand(RANDSEED);

    save_meta();
    prod_neighbors();

#ifdef RUN_SPEC_PARAM
    int param_idx;
    for (param_idx = 0; param_idx < sizeof(SPEC_PARAM_LIST)/(sizeof(double) * N_PARAM); param_idx ++) {
        b = SPEC_PARAM_LIST[param_idx][0];
        run(fp_aver);
    }
#else
    for (b = b_start; b < b_end || fabs(b - b_end) < 1e-7; b += b_interval)
        run(fp_aver);
#endif
    fclose(fp_aver);
    return 0;
}

