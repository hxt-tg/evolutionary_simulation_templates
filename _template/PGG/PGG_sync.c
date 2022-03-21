#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <error.h>
#include <math.h>
#include <time.h>

/* ====================== Settings ======================= */
  /* Global settings */
#define REC_SNAPSHOT  1             /* Not euqal 0 if record snapshots */
#define REC_STEPS     1             /* Not euqal 0 if record data per step */
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
#define TRY_TIMES     3             /* Average result in tried times */
  /* Running parameters */
double r_start = 3, r_end = 4, r_interval = 0.1;
double cost = 1;
/* =================== End of settings =================== */


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
    _fprint_int(fp, _rand, "Random seed");
    _fprint_int(fp, MC_STEPS, "Monte Carlo steps");
    _fprint_int(fp, AVG_STEPS, "Recorded MC steps at the end");
    _fprint_int(fp, TRY_TIMES, "Try times per parameter");
    fprintf(fp, "r=%g~%g with interval %g.\n", r_start, r_end, r_interval);
    fclose(fp);
}
/* ============== End of utility functions =============== */








double r;
int net[SIZE][N_NEIGH];   /* Players neighbor relations */
char stra[SIZE], last_stra[SIZE];
double last_payoff[SIZE];

/* Properties counter */
typedef enum {
    _C,
    _D,
    N_STRA
} Strategy;

int stra_cnt[N_STRA] = {0};

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
    for (i = 0; i < SIZE; i++) stra_cnt[stra[i] = randi(N_STRA)]++;
}

void calc_last_payoff_by_group() {
    double pool[SIZE];
    int x, i;
    for (x = 0; x < SIZE; x++) {
        pool[x] = stra[x] == _C ? cost/(N_NEIGH+1) : 0;
        for (i = 0; i < N_NEIGH; i++)
            pool[x] += stra[net[x][i]] == _C ? cost/(N_NEIGH+1) : 0;
        pool[x] *= r;
    }
    for (x = 0; x < SIZE; x++) {
        last_payoff[x] = pool[x]/(N_NEIGH+1);
        for (i = 0; i < N_NEIGH; i++)
            last_payoff[x] += pool[net[x][i]]/(N_NEIGH+1);
        last_payoff[x] -= stra[x] == _C ? cost : 0;
    }
}

void update_stra(int x) {
    int y = net[x][(int) randi(N_NEIGH)];
    /* Update strategy, if same then exit */
    if (last_stra[x] == last_stra[y]) return ;
    if (randf() < 1/(1+exp((last_payoff[x] - last_payoff[y])/K))){
        stra_cnt[stra[x]]--;
        stra[x] = last_stra[y];
        stra_cnt[stra[x]]++;
    }
}

/* Snapshot function temporarily not used */
void snapshot(int try_time, int step) {
    if (!REC_SNAPSHOT) return ;
    if (TRY_TIMES > 1) sprintf(_file_path, "output/snapshots/r=%.2lf_step=%d[%d].csv", r, step+1, try_time);
    else sprintf(_file_path, "output/snapshots/r=%.2lf_step=%d.csv", r, step+1);
    FILE *fp = fopen(_file_path, "w");
    int x;
    for (x = 0; x < SIZE; x++)
        fprintf(fp, "%d%c", stra[x], (x+1)%L?',':'\n');
    fclose(fp);
}

void run(FILE *fp_aver) {
    FILE *fp_steps = NULL;
    if (REC_STEPS) {
        sprintf(_file_path, "output/steps/r=%.2lf.csv", r);
        fp_steps = fopen(_file_path, "w");
        if (TRY_TIMES > 1) fprintf(fp_steps, "try,");
        fprintf(fp_steps, "step,f_C,f_D\n");
    }
    int i, x, step, rec, ta[N_STRA]={0}, try_time;      /* total amount for stras in all REC steps */
    double f[N_STRA]={0}, af[TRY_TIMES][N_STRA]={0};    /* freq and avg_freq */
    for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
        for (i = 0; i < N_STRA; i++) ta[i] = 0;
        rec = 0;
        init_randomly();
        snapshot(try_time, -1);
        for (step = 0; step < MC_STEPS; step++) {
            calc_last_payoff_by_group();
            for (x = 0; x < SIZE; x++)
                last_stra[x] = stra[x];
                
            for (x = 0; x < SIZE; x++)
                update_stra(x);
            
            for (i = 0; i < N_STRA; i++)
                f[i] = stra_cnt[i]/(double)SIZE;
                
            if (step > MC_STEPS-REC_STEPS-1) {
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
        printf(" r: %.2lf   avg: [%6.2lf%%,%6.2lf%%]       \n", r,
                af[try_time][_C]*100, af[try_time][_D]*100);
        
    }
    if (REC_STEPS) fclose(fp_steps);
    
    double aaf[N_STRA] = {0};       /* Average in all tries */
    for (i = 0; i < N_STRA; i++) {
        for (x = 0; x < TRY_TIMES; x ++)
            aaf[i] += af[x][i];
        aaf[i] /= TRY_TIMES;
    }
    fprintf(fp_aver, "%.2lf,%lf,%lf\n", r, aaf[_C], aaf[_D]);
    if (TRY_TIMES > 1) printf("[%d tries] r: %g  avg: [%6.2lf%%,%6.2lf%%]\n\n", 
            TRY_TIMES, r, aaf[_C]*100, aaf[_D]*100);
}

int main() {
    mkdir("output");
    if (REC_STEPS) mkdir("output/steps");
    if (REC_SNAPSHOT) mkdir("output/snapshots");
    
    /* Flush file */
    FILE *fp_aver = fopen("output/average.csv", "w");
    fprintf(fp_aver, "r,f_C,f_D\n");
    
    _rand = RANDSEED == -1 ? time(NULL) : RANDSEED;
    sgenrand(RANDSEED);
    
    save_meta();
    prod_neighbors();
    
    for (r = r_start; r < r_end || fabs(r - r_end) < 1e-7; r += r_interval)
        run(fp_aver);
    
    fclose(fp_aver);
    return 0;
}

