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


#include "cimnet/network.h"
#include "cimnet/random.h"


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

/* Running under special parameters. (Enabled by RUN_SPEC_PARAM) */
double SPEC_PARAM_LIST[][N_PARAM] = {
    /* {b}, */
    {1.00},
    {1.02},
    {1.06},
    {1.20},
};


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
    std::cerr << "Unsupport platform. (Used in make_dir)" << std::endl;
    exit(ENOTSUP);
#endif

    if (status && errno != EEXIST) {
        std::cerr << "Cannot mkdir 'output'." << std::endl;
        exit(errno);
    }
}
/* ============== End of utility functions =============== */



class PDGSimulation {
public:
    typedef enum {
        _C,
        _D,
        N_STRA
    } Strategy;
    
    typedef struct {
        Strategy stra;
    } NodeData;
    
    typedef None EdgeData;
    typedef Network<int, NodeData, None> PDGNetwork;
    /* Properties counter */
    int stra_cnt[N_STRA] = {0};
    /* Payoff matrix and its update (N_STRA*N_STRA matrix) */
    double payoff_matrix[N_STRA][N_STRA] = {
        {1, 0},
        {0, 0}
    };
    /* ST_MAT: Configure strategies matrix for initialization. */
    #define ST_MAT_ROW  2
    #define ST_MAT_COL  2
    const Strategy ST_MAT[ST_MAT_ROW][ST_MAT_COL] = {
        { _C, _D},
        { _D, _C},
    };
    
    /* FUNCTIONS */
    PDGSimulation(PDGNetwork &_net, double _b)
        : try_time(0), step(-1), b(_b), net(_net) {
        payoff_matrix[1][0] = b;
    }
    
    void init_randomly() {
        int i;
        for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
        for (i = 0; i < SIZE; i++) stra_cnt[net[i].stra = (Strategy)randi(N_STRA)]++;
    }

    void init_by_st_mat() {
    #ifndef GRID_NET
        sprintf(stderr, "Matrix initialization is only supported on GridNetwork.\n");
        exit(1);
    #endif
        int i, j, block_row = L/ST_MAT_ROW, block_col = L/ST_MAT_COL;
        for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
        for (i = 0; i < L; i++)
            for (j = 0; j < L; j++) 
                stra_cnt[net[i*L+j].stra = ST_MAT[i/block_row][j/block_col]]++;
    }

    double payoff(int x) {
        /* Isolated node */
        if (net.degree(x) == 0) return 0;
        double pay = 0;
        for (auto &n : net.neighbors(x))
            pay += payoff_matrix[net[x].stra][net[n].stra];
        return pay;
    }
    
    void update_stra(int x) {
        /* Isolated node */
        if (net.degree(x) == 0) return ;
        int y = net.random_neighbor(x);
        /* Update strategy, if same then exit */
        if (net[x].stra == net[y].stra) return ;
        if (randf() < 1/(1+exp((payoff(x) - payoff(y))/K))){
            stra_cnt[net[x].stra]--;
            net[x].stra = net[y].stra;
            stra_cnt[net[y].stra]++;
        }
    }
    
    /* Snapshot function temporarily not used */
    void matrix_snapshot() {
    #ifndef REC_SNAPSHOT
        return ;
    #endif
    #ifndef GRID_NET
        sprintf(ftderr, "Matrix snapshot is only supported on GridNetwork.\n");
        exit(1);
    #endif
        if (TRY_TIMES > 1) sprintf(_file_path, "output/snapshots/b=%.2lf_step=%d[%d].csv", b, step+1, try_time);
        else sprintf(_file_path, "output/snapshots/b=%.2lf_step=%d.csv", b, step+1);
        FILE *fp = fopen(_file_path, "w");
        int x;
        for (x = 0; x < SIZE; x++)
            fprintf(fp, "%d%c", net[x].stra, (x+1)%L?',':'\n');
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
        int i, x, rec, ta[N_STRA]={0};      /* total amount for stras in all REC steps */
        double f[N_STRA]={0}, af[TRY_TIMES][N_STRA]={0};    /* freq and avg_freq */
        for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
            for (i = 0; i < N_STRA; i++) ta[i] = 0;
            rec = 0;
    #ifdef USE_MATRIX_INIT
            init_by_st_mat();
    #else
            init_randomly();
    #endif
            matrix_snapshot();
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
    
private:
    int try_time;
    int step;
    double b;
    PDGNetwork net;
};




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


    std::cout << "Construct network ..." << std::endl;
#if defined(GRID_NET)
    GridNetwork<PDGSimulation::NodeData, PDGSimulation::EdgeData> net(L, L, N_NEIGH);
#else
    std::cerr << "Unspecified net structure." << std::endl;
    exit(ENOTSUP);
#endif
    std::cout << net << std::endl;

#ifdef RUN_SPEC_PARAM
    int param_idx;
    for (param_idx = 0; param_idx < sizeof(SPEC_PARAM_LIST)/(sizeof(double) * N_PARAM); param_idx ++) {
        b = SPEC_PARAM_LIST[param_idx][0];
        PDGSimulation pdg(net, b);
        pdg.run(fp_aver);
    }
#else
    for (double b = b_start; b < b_end || fabs(b - b_end) < 1e-7; b += b_interval) {
        PDGSimulation pdg(net, b);
        pdg.run(fp_aver);
    }
#endif
    fclose(fp_aver);
    return 0;
}

