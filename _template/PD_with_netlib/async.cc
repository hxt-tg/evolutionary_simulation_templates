#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#if defined(__linux__)
#include <sys/stat.h>
#include <linux/limits.h>
#elif defined(__MINGW32__) || defined(__MINGW64__)
#include <sys/stat.h>
#include <climits>
#else              // defined(_WIN32) || defined(_WIN64)
#include <climits>
#endif


#include "netlib/common_net.h"
#include "netlib/sfnet.h"
#include "netlib/random.h"


/* ====================== Settings ======================= */
/* Global settings */
#define REC_SNAPSHOT  1             /* Not euqal 0 if record snapshots */
#define REC_STEPS     1             /* Not euqal 0 if record data per step */
#define RANDSEED      3145215       /* Fixed random seed, use time if -1 */
/* Network settings */
#define GRID_NET

#ifdef GRID_NET
#define L             100           /* Side length of a network */
#define SIZE          (L * L)       /* Total number of nodes */
#define N_NEIGH       4             /* Number of neighbors for every node */
#endif
/* Monte Carlo settings */
#define K             0.1           /* Fermi coefficient */
#define MC_STEPS      10000         /* Number of MC steps */
#define AVG_STEPS     5000          /* Number of steps at the end for averaging */
#define REFRESH_FREQ  100           /* Control the frequency of refresh screen, 1 for all the time */
#define TRY_TIMES     1             /* Average result in tried times */
/* Running parameters */
double b_start = 1, b_end = 1.5, b_interval = 0.05;
/* =================== End of settings =================== */

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
    fprintf(fp, "b=%g~%g with interval %g.\n", b_start, b_end, b_interval);
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
UndirectedNet *net;

/* Properties counter */
typedef enum {
    _C,
    _D,
    N_STRA
} Strategy;

Strategy stra[SIZE];
int stra_cnt[N_STRA] = {0};

/* Payoff matrix and its update (N_STRA*N_STRA matrix) */
double payoff_matrix[N_STRA][N_STRA] = {{1, 0},
    {0, 0}};
/* Call update_matrix(b) after loop for b */
#define update_matrix(b) do {    \
    payoff_matrix[1][0]=b;       \
} while (0)

void build_net() {
#if defined(GRID_NET)
    net = new GridNet(L, L, N_NEIGH);
#else
    std::cerr << "Unspecified net structure." << std::endl;
    exit(ENOTSUP);
#endif
    net->buildNeighbors();
}

void init_randomly() {
    int i;
    for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
    for (i = 0; i < SIZE; i++) stra_cnt[stra[i] = (Strategy)randi(N_STRA)]++;
}


double payoff(int x) {
    /* Isolated node */
    if (net->vexDegree(x) == 0) return 0;
    double pay = 0;
    NeiArr nei = net->getNeighbors(x);
    int i;
    for (i = 0; i < net->vexDegree(x); i++)
        pay += payoff_matrix[stra[x]][stra[nei[i]]];
    return pay;
}

void update_stra(int x) {
    /* Isolated node */
    if (net->vexDegree(x) == 0) return;
    NeiArr nei = net->getNeighbors(x);
    int y = nei[(int) randi (net->vexDegree(x))];
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
    if (!REC_SNAPSHOT) return ;
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
    if (REC_STEPS) {
        sprintf(_file_path, "output/steps/b=%.2lf.csv", b);
        fp_steps = fopen(_file_path, "w");
        if (TRY_TIMES > 1) fprintf(fp_steps, "try,");
        fprintf(fp_steps, "step,f_C,f_D\n");
    }
    int i, x, step, rec, ta[N_STRA]={0}, try_time;      /* total amount for stras in all REC steps */
    double f[N_STRA]={0}, af[TRY_TIMES][N_STRA]={0};    /* freq and avg_freq */
    update_matrix(b);
    for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
        for (i = 0; i < N_STRA; i++) ta[i] = 0;
        rec = 0;
        init_randomly();
        snapshot(try_time, -1);
        for (step = 0; step < MC_STEPS; step++) {
            for (x = 0; x < SIZE; x++)
                update_stra((int)randi(SIZE));

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
        printf(" b: %.2lf   avg: [%6.2lf%%,%6.2lf%%]       \n", b,
                af[try_time][_C]*100, af[try_time][_D]*100);

    }
    if (REC_STEPS) fclose(fp_steps);

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
    if (REC_STEPS) make_dir("output/steps");
    if (REC_SNAPSHOT) make_dir("output/snapshots");

    /* Flush file */
    FILE *fp_aver = fopen("output/average.csv", "w");
    fprintf(fp_aver, "b,f_C,f_D\n");

    _rand = RANDSEED == -1 ? time(NULL) : RANDSEED;
    sgenrand(RANDSEED);

    save_meta();
    build_net();

    for (b = b_start; b < b_end || fabs(b - b_end) < 1e-7; b += b_interval)
        run(fp_aver);

    fclose(fp_aver);
    return 0;
}

