#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <errno.h>     /* for E* */
#include <unistd.h>    /* for mkdir */
#include <getopt.h>    /* for param read */
#include <sys/types.h> /* for S_* */


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

/* Network settings */
#define GRID_NET   /* Temporarily GRID_NET only */

/* Constants */
#define NEINUM                4
#define NAME_MAX            255
char *progname=NULL;

/* Buffers */
char _file_path[PATH_MAX];

/* Optional parameters */
#define DFT_K                 0.1
#define DFT_L                 100
#define DFT_MC_STEPS        10000
#define DFT_REC_STEPS        1000
#define DFT_TRY_TIMES           3
#define DFT_REFRESH_FREQ      100
#define DFT_SEED          3145215

double K=DFT_K;
bool NO_FILE_OUTPUT=false,
     NO_SNAPSHOT=false;
int L=DFT_L,
    SIZE=DFT_L*DFT_L,
    SEED=DFT_SEED,
    MC_STEPS=DFT_MC_STEPS,
    REC_STEPS=DFT_REC_STEPS,
    TRY_TIMES=DFT_TRY_TIMES,
    REFRESH_FREQ=DFT_REFRESH_FREQ;

/* Global path name */
char OUTPUT_DIR[NAME_MAX]="",
     PARAM_RANGE_NAME[NAME_MAX]="",
     META_FILE[NAME_MAX]="",
     AVG_FILE[NAME_MAX]="",
     STEPS_DIR[NAME_MAX]="",
     SNAPSHOTS_DIR[NAME_MAX]="";

/* Required parameters */
double b_start = 1, b_end = 1.2, b_interval = 0.05;

const struct option long_options[] = {
    {"temptation",   required_argument, 0, 'b' },
    {"side-length",  required_argument, 0, 'L' },
    {"temperature",  required_argument, 0, 'K' },
    {"refresh-freq", required_argument, 0, 'F' },
    {"seed",         required_argument, 0, 0x1001 },
    {"try-times",    required_argument, 0, 0x1002 },
    {"mc-steps",     required_argument, 0, 0x1003 },
    {"record-steps", required_argument, 0, 0x1004 },
    {"no-file",      no_argument,       0, 0x1005 },
    {"no-snapshot",  no_argument,       0, 0x1006 },
    {"help",         no_argument,       0, 'h'    },
    {0,              0,                 0,  0     }
};

/* Parameters processor */
const char *_fmtss        = "    -%c  --%-14s%-32s\"%s\"\n",      /* short string */
           *_fmtsi        = "    -%c  --%-14s%-32s %d\n",         /* short int */
           *_fmtsd        = "    -%c  --%-14s%-32s %g\n",         /* short double */
           *_fmtsdr       = "    -%c  --%-14s%-32s <Required>\n", /* short double (required) */
           *_fmtls        = "        --%-14s%-32s\"%s\"\n",       /* long string */
           *_fmtli        = "        --%-14s%-32s %d\n",          /* long int */
           *_fmtln        = "        --%-14s%-32s\n",             /* long (no-argument) */
           *_fmt_desc     = "                          %-30s\n";
void show_param_help() {
    static int is_showed = 0;
    if (is_showed) return ;
    is_showed++;
    printf(
"Usage: .\\%s -b b_start,b_end,b_step [optinal params]\n"
"  Each optional parameter has a short and long name.\n"
"    (e.g. \"-b\" or \"--temptation\")\n"
"  They have the same meaning. Use one of them if set.\n"
"[OPTIONS]\n", progname);
    printf(" short  long            Param description               Default value\n");
    printf(_fmtsdr,  'b', "temptation",   "PD temptation value");
    printf(_fmtsi,   'L', "side-length",  "Side length of square lattice",  DFT_L);
    printf(_fmtsd,   'K', "temperature",  "Temperature coefficient",        DFT_K);
    printf(_fmtsi,   'F', "refresh-freq", "Steps of verbose refresh",       DFT_REFRESH_FREQ);
    printf(_fmt_desc,                     "at every \"refresh-freq\" steps");
    printf(_fmt_desc,                     "(if 0, turn off refresh)");
    printf(_fmtli,        "seed",         "Global random seed",             DFT_SEED);
    printf(_fmt_desc,                     "(This option is not used yet.)");
    printf(_fmtli,        "try-times",    "Try times for whole MCS",        DFT_TRY_TIMES);
    printf(_fmtli,        "mc-steps",     "Monte Carlo steps",              DFT_MC_STEPS);
    printf(_fmtli,        "record-steps", "Number of last rounds for avg",  DFT_REC_STEPS);
    printf(_fmtln,        "no-file",      "Turn off file output");
    printf(_fmtln,        "no-snapshot",  "Turn off network snapshot");
    printf(_fmtln,        "help",         "Show this help");
    printf(
"[NOTICE]\n"
"   There should not be any space in \"temptation\" option.\n"
"   If steps and average output are all disabled, output\n"
"   folder will not be created, so as metadata file.\n");
printf("[EXAMPLE] .\\%s -b1,1.2,0.05 --no-file\n", progname);
}

void die(int ecode, const char *einfo) {
    if (einfo)
        fprintf(stderr, "[ERR] %s\n", einfo);
    fprintf(stderr, "%s: Use \"--help\" option for more information\n", progname);
    exit(ecode);
}

#define _check_not_0(x) do {if (x < 1) die(EINVAL, "#x should be larger than 0.");} while(0)
void read_parameters(int argc, char **argv) {
    if (argc < 2) {
        show_param_help();
        fprintf(stderr, "\n[ERR] This program requires parameters.\n"
                        "  You may need run it in command line. (or terminal)\n");
        exit(EINVAL);
    }
    int c, has_b=0;
    while (1) {
        int option_index = 0;
        c = getopt_long(argc, argv, "b:L:K:F:h", long_options, &option_index);
        if (c == -1) break;
        switch (c) {
        case 'b':
            if (sscanf(optarg, "%lf,%lf,%lf", &b_start, &b_end, &b_interval) < 1)
                die(EINVAL, "Invalid temptation (b).");
            has_b = 1;
            break;
        case 'L':
            if (sscanf(optarg, "%d", &L) < 1)
                die(EINVAL, "Invalid side length of square lattice.");
            break;
        case 'K':
            if (sscanf(optarg, "%lf", &K) < 1)
                die(EINVAL, "Invalid temperation coefficient.");
            break;
        case 'F':
            if (sscanf(optarg, "%d", &REFRESH_FREQ) < 1)
                die(EINVAL, "Invalid refresh frequency.");
            break;
        case 0x1001:
            if (sscanf(optarg, "%d", &SEED) < 1)
                die(EINVAL, "Invalid random seed.");
            break;
        case 0x1002:
            if (sscanf(optarg, "%d", &TRY_TIMES) < 1)
                die(EINVAL, "Invalid try times.");
            break;
        case 0x1003:
            if (sscanf(optarg, "%d", &MC_STEPS) < 1)
                die(EINVAL, "Invalid Monte Carlo steps.");
            break;
        case 0x1004:
            if (sscanf(optarg, "%d", &REC_STEPS) < 1)
                die(EINVAL, "Invalid record steps.");
            break;
        case 0x1005:
            NO_FILE_OUTPUT = true;
            break;
        case 0x1006:
            NO_SNAPSHOT = true;
            break;
        case 'h':
            show_param_help();
            exit(0);
        case '?':
            fprintf(stderr, "%s: Use \"--help\" option for more information\n", progname);
            exit(EINVAL);
        default:
            die(EINVAL, "Unknown parameter error.");
        }
    }
    
    if (!has_b) die(EINVAL, "Option temptation (b) is nessesary.");
    if (optind < argc) {
       fprintf(stderr, "[WARN] non-options: ");
       while (optind < argc)
           fprintf(stderr, "%s ", argv[optind++]);
       fprintf(stderr, "\n");
    }
    /* Check vars */
    if (L < 2) die(EINVAL, "L should be larger than 1.");
    if (K < 0 || K > 1) die(EINVAL, "K should be in range [0, 1].");
    if (TRY_TIMES < 1) die(EINVAL, "TRY_TIMES should be larger than 0.");
    if (MC_STEPS < 1) die(EINVAL, "MC_STEPS should be larger than 0.");
    if (REC_STEPS < 1) die(EINVAL, "REC_STEPS should be larger than 0.");
    if (REC_STEPS > MC_STEPS) die(EINVAL, "REC_STEPS should be larger than MC_STEPS.");
    if (b_start > b_end) die(EINVAL, "b_start should be lower than b_end.");
    if (fabs(b_interval) < 1e-7 || b_interval < 0) die(EINVAL, "b_interval should be larger than 0.");
    /* Construct initial values of global vars */
    SIZE = L*L;
    if (SEED != -1) sgenrand(SEED);
    else sgenrand((unsigned int)time(NULL));
    sprintf(PARAM_RANGE_NAME, "b=%g~%g", b_start, b_end);
    if (sprintf(OUTPUT_DIR, "output_%s", PARAM_RANGE_NAME) > NAME_MAX-1)
        die(ENAMETOOLONG, "Output directory name is too long. (max=254)");
    if (sprintf(STEPS_DIR, "%s/steps", OUTPUT_DIR) > NAME_MAX-1)
        die(ENAMETOOLONG, "Steps directory name is too long. (max=254)");
    if (sprintf(SNAPSHOTS_DIR, "%s/snapshots", OUTPUT_DIR) > NAME_MAX-1)
        die(ENAMETOOLONG, "Snapshots directory name is too long. (max=254)");
    if (sprintf(META_FILE, "%s/_meta.csv", OUTPUT_DIR) > NAME_MAX-1)
        die(ENAMETOOLONG, "Meta file name is too long. (max=254)");
    if (sprintf(AVG_FILE, "%s/average_%s.csv", OUTPUT_DIR, PARAM_RANGE_NAME) > NAME_MAX-1)
        die(ENAMETOOLONG, "Average file name is too long. (max=254)");
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
        fprintf(stderr, "Cannot mkdir '%s'.\n", path);
        exit(errno);
    }
}

#define _fprint_str(fp, x) fprintf(fp, "%s,\"%s\"\n", #x, x)
#define _fprint_int(fp, x) fprintf(fp, "%s,%d\n", #x, x)
#define _fprint_dbl(fp, x) fprintf(fp, "%s,%g\n", #x, x)

void save_metadata() {
    if (NO_FILE_OUTPUT) return ;
    make_dir(OUTPUT_DIR);
    FILE *fp = fopen(META_FILE, "w");
    _fprint_dbl(fp, b_start);
    _fprint_dbl(fp, b_end);
    _fprint_dbl(fp, b_interval);
    _fprint_int(fp, L);
    _fprint_int(fp, SIZE);
    _fprint_dbl(fp, K);
    _fprint_int(fp, SEED);
    _fprint_int(fp, MC_STEPS);
    _fprint_int(fp, REC_STEPS);
    _fprint_int(fp, TRY_TIMES);
    _fprint_str(fp, AVG_FILE);
    _fprint_str(fp, STEPS_DIR);
    fclose(fp);
}


double b;

/* Properties counter */
typedef enum {
    _C,
    _D,
    N_STRA
} Strategy;

const char *STRATEGY_TAGS[N_STRA] = {
    "Cooperate",
    "Defect"
};

typedef struct {
    Strategy stra;
} NodeData;

typedef None EdgeData;

typedef Network<int, NodeData, EdgeData> PDNetwork;

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

void init_randomly(PDNetwork &net) {
    int i;
    for (i = 0; i < N_STRA; i++) stra_cnt[i] = 0;
    for (i = 0; i < SIZE; i++) stra_cnt[net[i].stra = (Strategy)randi(N_STRA)]++;
}

double payoff(PDNetwork &net, int x) {
    /* Isolated node */
    if (net.degree(x) == 0) return 0;
    double pay = 0;
    for (auto &n : net.neighbors(x))
        pay += payoff_matrix[net[x].stra][net[n].stra];
    return pay;
}

void update_stra(PDNetwork &net, int x) {
    /* Isolated node */
    if (net.degree(x) == 0) return ;
    auto x_nei = net.neighbors(x);
    int y = x_nei[randi(x_nei.size())];
    /* Update strategy, if same then exit */
    if (net[x].stra == net[y].stra) return ;
    if (randf() < 1/(1+exp((payoff(net, x) - payoff(net, y))/K))){
        stra_cnt[net[x].stra]--;
        net[x].stra = net[y].stra;
        stra_cnt[net[y].stra]++;
    }
}

/* Snapshot (grid network) */
void snapshot(PDNetwork &net, int try_time, int step) {
    if (NO_FILE_OUTPUT || NO_SNAPSHOT) return ;
#ifndef GRID_NET
    fprintf(stderr, "Snapshot is only supported on GridNetwork.\n");
    exit(-1);
#endif
    make_dir(SNAPSHOTS_DIR);
    if (TRY_TIMES > 1) sprintf(_file_path, "%s/b=%.2lf_step=%d[%d].csv", 
        SNAPSHOTS_DIR, b, step+1, try_time);
    else sprintf(_file_path, "%s/b=%.2lf_step=%d.csv", 
        SNAPSHOTS_DIR, b, step+1);
    FILE *fp = fopen(_file_path, "w");
    int x;
    for (x = 0; x < SIZE; x++)
        fprintf(fp, "%c%c", STRATEGY_TAGS[net[x].stra][0], (x+1)%L?',':'\n');
    fclose(fp);
}

void run_param(PDNetwork &net) {
    update_matrix(b);
    int step, rec, ta[N_STRA], try_time;          /* total amount for stras in all REC steps */
    double f[N_STRA]={0}, af[TRY_TIMES][N_STRA];  /* freq and avg_freq */
    FILE *fp_avg = NULL, *fp_steps = NULL;
    if (!NO_FILE_OUTPUT) {
        fp_avg = fopen(AVG_FILE, "a");
        sprintf(_file_path, "%s/b=%.2lf.csv", STEPS_DIR, b);
        fp_steps = fopen(_file_path, "w");
        fprintf(fp_steps, "%sstep", (TRY_TIMES > 1 ? "try_time," : ""));
        for (int i = 0; i < N_STRA; i++)
            fprintf(fp_steps, ",f_%c", STRATEGY_TAGS[i][0]);
        fprintf(fp_steps, "\n");
    }
    for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
        for (int i = 0; i < N_STRA; i++) af[try_time][i] = ta[i] = 0;
        rec = 0;

        init_randomly(net);
        snapshot(net, try_time, -1);

        for (step = 0; step < MC_STEPS; step++) {
            for (int x = 0; x < SIZE; x ++)
                update_stra(net, (int)randi(SIZE));
            
            // Strategies stats
            for (int i = 0; i < N_STRA; i++)
                f[i] = stra_cnt[i]/(double)SIZE;
            if (step > MC_STEPS-REC_STEPS-1) {
                for (int i = 0; i < N_STRA; i++)
                    ta[i] += stra_cnt[i];
                rec++;
            }
            
            // Write steps file
            if (!NO_FILE_OUTPUT) {
                if (TRY_TIMES > 1) fprintf(fp_steps, "%d,", try_time);
                fprintf(fp_steps, "%d", step);
                for (int i = 0; i < N_STRA; i++)
                    fprintf(fp_steps, ",%g", f[i]);
                fprintf(fp_steps, "\n");
            }
            
            // Update screen
            if (REFRESH_FREQ && step % REFRESH_FREQ == 0) {
                if (TRY_TIMES > 1) printf("\r [%d]", try_time);
                else printf("\r");
                printf(" Step: %6d  f: [", step);
                for (int i = 0; i < N_STRA; i++)
                    printf("%s%6.2lf%%", (i == 0 ? "": ","), f[i] * 100);
                printf("]");
            }
            
            // Pre-break when not evolving
            if (stra_cnt[_C] == SIZE || stra_cnt[_D] == SIZE) {
                if (step++ < MC_STEPS-REC_STEPS)
                    for (int i = 0; i < N_STRA; i++)
                        af[try_time][i] = stra_cnt[i] ? 1 : 0;
                break;
            };
        }
        // Average in one try
        if (rec)
            for (int i = 0; i < N_STRA; i++)
                af[try_time][i] = (double)ta[i]/(SIZE * rec);
        
        // Update average stats
        if (TRY_TIMES > 1) printf("\r [%d]", try_time);
        else printf("\r");
        printf(" b: %.2f   avg: [", b);
        for (int i = 0; i < N_STRA; i++)
            printf("%s%6.2lf%%", (i == 0 ? "": ","), af[try_time][i]*100);
        printf("]        \n");
    }
    // Average in all tries
    double aaf[N_STRA] = {0};
    for (int i = 0; i < N_STRA; i++) {
        for (int t = 0; t < TRY_TIMES; t++)
            aaf[i] += af[t][i];
        aaf[i] /= TRY_TIMES;
    }

    // Write average file and update to stdout
    fprintf(fp_avg, "%.2lf", b);
    for (int i = 0; i < N_STRA; i++)
        fprintf(fp_avg, ",%lf", aaf[i]);
    fprintf(fp_avg, "\n");
    if (TRY_TIMES > 1) {
        printf("[%d tries] b: %.2f   avg: [",TRY_TIMES, b);
        for (int i = 0; i < N_STRA; i++)
            printf("%s%6.2lf%%", (i == 0 ? "": ","), aaf[i]*100);
        printf("]\n\n");
    }

    if (!NO_FILE_OUTPUT) {
        fclose(fp_steps);
        fclose(fp_avg);
    }
}

char *get_progname(char ** const argv) {
    char *slash_pos = strrchr(argv[0], '/');
    if (!slash_pos) slash_pos = strrchr(argv[0], '\\');
    return slash_pos ? slash_pos+1 : argv[0];
}

int main(int argc, char **argv) {
    progname = get_progname(argv);
    read_parameters(argc, argv);
    save_metadata();
    
    /* Flush average file */
    if (!NO_FILE_OUTPUT) {
        make_dir(OUTPUT_DIR);
        make_dir(STEPS_DIR);
        FILE *fp_avg = fopen(AVG_FILE, "w");
        fprintf(fp_avg, "b");
        for (int i = 0; i < N_STRA; i++)
            fprintf(fp_avg, ",f_%c", STRATEGY_TAGS[i][0]);
        fprintf(fp_avg, "\n");
        fclose(fp_avg);
    }

    std::cout << "Construct grid network ..." << std::endl;
    GridNetwork<NodeData, EdgeData> net(L, L, NEINUM);
    std::cout << net << std::endl;

	for (b = b_start; b < b_end || fabs(b - b_end) < 1e-7; b += b_interval)
        run_param(net);
	
    return 0;
}

