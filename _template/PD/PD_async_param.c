#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <errno.h>     /* for E* */
#include <unistd.h>    /* for mkdir */
#include <getopt.h>    /* for param read */
#include <sys/types.h> /* for S_* */
#include <string.h>
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

/* Random function (Do not modify below) */
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
/* End of random (Do not modify above) */




/* Constants */
#define NEINUM                4
#define STRANUM               2
#define NAME_MAX            255
#define META_FN          "_meta"
#define SNAPSHOT_DIR "snapshots"
char *progname="";

/* Buffers */
char _file_path[PATH_MAX];

/* Optional parameters */
#define DFT_K                 0.1
#define DFT_L                 100
#define DFT_MC_STEPS        10000
#define DFT_TRY_TIMES           3
#define DFT_REC_STEPS        1000
#define DFT_REFRESH_FREQ      100            /* Control the frequency of refresh screen, 1 for all the time */
#define DFT_OUTPUT_DIR    "output"
#define DFT_STEPS_DIR      "steps"
#define DFT_AVG_FN       "average"
#define DFT_STRA_STR          "CD"

double K=DFT_K;

int L=DFT_L,
    SIZE=DFT_L*DFT_L,
    SEED=DFT_SEED,
    MC_STEPS=DFT_MC_STEPS,
    REC_STEPS=DFT_REC_STEPS,
    TRY_TIMES=DFT_TRY_TIMES,
    REFRESH_FREQ=DFT_REFRESH_FREQ;

char OUTPUT_DIR[NAME_MAX]=DFT_OUTPUT_DIR,
     STEPS_DIR[NAME_MAX]=DFT_STEPS_DIR,
     AVG_FN[NAME_MAX]=DFT_AVG_FN,
     STRA_STR[STRANUM+1]=DFT_STRA_STR;

/* Required parameters */
double bs, be, bi;        /* b, start, end, interval */

const struct option long_options[] = {
    {"output-dir",   required_argument, 0,  0 },
    {"steps-dir",    required_argument, 0,  0 },
    {"avg-file",     required_argument, 0,  0 },
    {"random-seed",  required_argument, 0,  0 },
    {"stra-name",    required_argument, 0,  0 },
    {"side-length",  required_argument, 0, 'L'},
    {"temperature",  required_argument, 0, 'K'},
    {"try-times",    required_argument, 0, 'T'},
    {"mc-steps",     required_argument, 0, 'M'},
    {"record-steps", required_argument, 0, 'R'},
    {"refresh-freq", required_argument, 0,  0 },
    {"temptation",   required_argument, 0, 'b'},
    {"help",         no_argument,       0,  0 },
    {0,              0,                 0,  0 }
};

/* Parameters processor */
const char *_fmtss        = "    -%c  --%-14s%-32s\"%s\"\n",
           *_fmtsi        = "    -%c  --%-14s%-32s %d\n",
           *_fmtsd        = "    -%c  --%-14s%-32s %g\n",
           *_fmtls        = "        --%-14s%-32s\"%s\"\n",
           *_fmtli        = "        --%-14s%-32s %d\n",
           *_fmt_desc     = "                          %-30s\n";
void show_param_help() {
    static int is_showed = 0;
    if (is_showed) return ;
    is_showed++;
    printf(
"Usage: %s -b b_start,b_end,b_step [optinal params]\n"
"  Each optional parameter has a short and long name.\n"
"    (e.g. \"-b\" or \"--temptation\")\n"
"  They have the same meaning. Use one of them if set.\n"
" [OPTIONS]\n", progname);
    printf(" short  long            Param description               Default value\n");
    printf(_fmtls, "output-dir",   "Output folder name",            DFT_OUTPUT_DIR);
    printf(_fmtls, "steps-dir",    "Steps output folder name",      DFT_STEPS_DIR);
    printf(_fmt_desc,  "(if \"*\", no steps output)");
    printf(_fmtls, "avg-file",     "Average output file name",      DFT_AVG_FN);
    printf(_fmt_desc,  "(if \"*\", no average output)");
    printf(_fmtls, "stra-name",    "Strategies name tag",           DFT_STRA_STR);
    printf(_fmt_desc,  "(This option is not used yet.)");
    printf(_fmtli, "random-seed",  "Global random seed",            DFT_SEED);
    printf(_fmt_desc,  "(if -1, use sys time as seed)");
    printf(_fmtli, "refresh-freq", "Steps of verbose refresh",      DFT_REFRESH_FREQ);
    printf(_fmt_desc,  "at every \"refresh-freq\" steps");
    printf(_fmt_desc,  "(if 0, turn off refresh)");
    printf(_fmtsi, 'L', "side-length",  "Side length of square lattice", DFT_L);
    printf(_fmtsd, 'K', "temperature",  "Temperature coefficient",       DFT_K);
    printf(_fmtsi, 'T', "try-times",    "Try times for whole MCS",       DFT_TRY_TIMES);
    printf(_fmtsi, 'M', "mc-steps",     "Monte Carlo steps",             DFT_MC_STEPS);
    printf(_fmtsi, 'R', "record-steps", "Number of last rounds for avg", DFT_REC_STEPS);
    printf("        --%-14s%-32s\n", "help", "Show this help");
    printf(
" [NOTICE ]\n"
"   There should not be any space in \"temptation\" option.\n"
"   If steps and average output are all disabled, output\n"
"   folder will not be created, so as metadata file.\n");
printf(" [EXAMPLE]   .\\%s -b1,1.2,0.05 -s* -a* (No file output)\n", progname);
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
        c = getopt_long(argc, argv, "L:K:T:M:R:b:", long_options, &option_index);
        if (c == -1) break;
        switch (c) {
        case 'd':
            if (sprintf(OUTPUT_DIR, "%s", optarg) > NAME_MAX-1)
                die(ENAMETOOLONG, "Output folder name is too long. (max=254)");
            break;
        case 's':
            if (sprintf(STEPS_DIR, "%s", optarg) > NAME_MAX-1)
                die(ENAMETOOLONG, "Steps folder name is too long. (max=254)");
            if (strchr(STEPS_DIR, '*')) STEPS_DIR[0] = 0;
            break;
        case 'a':
            if (sprintf(AVG_FN, "%s", optarg) > NAME_MAX-1)
                die(ENAMETOOLONG, "Average file name is too long. (max=254)");
            if (strchr(AVG_FN, '*')) AVG_FN[0] = 0;
            break;
        case 'p':
            if (sprintf(VALUE_SEP, "%s", optarg) > SEP_MAX-1)
                die(EINVAL, 
                    "Separator should be less than 4. (One comma is recommanded)");
            break;
        case 'e':
            if (sprintf(DATA_FN_EXT, "%s", optarg) > EXT_MAX-1)
                die(EINVAL, 
                    "Data file path should be less than 32. (csv is recommanded)");
            break;
        case 'r':
            if (sscanf(optarg, "%d", &SEED) < 1)
                die(EINVAL, "Invalid random seed.");
            break;
        case 'S':
            if (sprintf(STRA_STR, "%s", optarg) != STRANUM)
                die(ENAMETOOLONG, "STRA_STR should contain exactly 2 chars.");
            break;
        case 'L':
            if (sscanf(optarg, "%d", &L) < 1)
                die(EINVAL, "Invalid side length of square lattice.");
            break;
        case 'K':
            if (sscanf(optarg, "%lf", &K) < 1)
                die(EINVAL, "Invalid temperation coefficient.");
            break;
        case 'T':
            if (sscanf(optarg, "%d", &TRY_TIMES) < 1)
                die(EINVAL, "Invalid try times.");
            break;
        case 'M':
            if (sscanf(optarg, "%d", &MC_STEPS) < 1)
                die(EINVAL, "Invalid Monte Carlo steps.");
            break;
        case 'R':
            if (sscanf(optarg, "%d", &REC_STEPS) < 1)
                die(EINVAL, "Invalid record steps.");
            break;
        case 'F':
            if (sscanf(optarg, "%d", &REFRESH_FREQ) < 1)
                die(EINVAL, "Invalid refresh frequency.");
            break;
        case 'b':
            if (sscanf(optarg, "%lf,%lf,%lf", &bs, &be, &bi) < 1)
                die(EINVAL, "Invalid temptation (b).");
            has_b = 1;
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
    if (L < 2) die(EINVAL, "L should be larger than 1.");
    if (K < 0 || K > 1) die(EINVAL, "K should be in range [0, 1].");
    if (TRY_TIMES < 1) die(EINVAL, "TRY_TIMES should be larger than 0.");
    if (MC_STEPS < 1) die(EINVAL, "MC_STEPS should be larger than 0.");
    if (REC_STEPS < 1) die(EINVAL, "REC_STEPS should be larger than 0.");
    if (REC_STEPS > MC_STEPS) die(EINVAL, "REC_STEPS should be larger than MC_STEPS.");
    if (bs > be) die(EINVAL, "b_start should be lower than b_end.");
    if (fabs(bi) < 1e-7 || bi < 0) die(EINVAL, "b_increse should be larger than 0.");
    if (strlen(STRA_STR) != STRANUM) die(EINVAL, "STRA_STR should contain 2 chars.");
    if (strchr(STRA_STR, '\n') || strchr(STRA_STR, '\r'))
        die(EINVAL, "STRA_STR should not contain '\\n' or '\\r'.");
    SIZE = L*L;
}

#define _fprint_str(fp, x) fprintf(fp, "%s,\"%s\"\n", #x, x)
#define _fprint_int(fp, x) fprintf(fp, "%s,%d\n", #x, x)
#define _fprint_dbl(fp, x) fprintf(fp, "%s,%g\n", #x, x)

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

void save_metadata() {
    sprintf(_file_path, "./%s", OUTPUT_DIR);
    make_dir(_file_path);
    sprintf(_file_path, "%s/%s.%s", OUTPUT_DIR, META_FN, DATA_FN_EXT);
    FILE *fp = fopen(_file_path, "w");
    _fprint_str(fp, VALUE_SEP);
    _fprint_str(fp, STRA_STR);
    _fprint_int(fp, L);
    _fprint_dbl(fp, K);
    _fprint_int(fp, SIZE);
    _fprint_int(fp, SEED);
    _fprint_int(fp, MC_STEPS);
    _fprint_int(fp, REC_STEPS);
    _fprint_int(fp, TRY_TIMES);
    fprintf(fp, "b_start,%g\nb_end,%g\nb_increase,%g\n", bs, be, bi);
    fprintf(fp, "STEPS_FILE,\"%s/%s.%s\"\nAVG_FILE,\"%s/%s.%s\"\n", 
        OUTPUT_DIR, AVG_FN, DATA_FN_EXT, OUTPUT_DIR, STEPS_DIR, DATA_FN_EXT);
    fclose(fp);
}

/* Global vars */
typedef int _neigh[NEINUM];
_neigh *net;       /* Players neighbor relations */

/* Properties counter */
char *stra;
int stra_cnt[STRANUM];    /* 0 for C and 1 for D */
double b;

/* Payoff matrix and its update */
double payoff_matrix[STRANUM][STRANUM] = {
    {1, 0},
    {0, 0}};
/* Call update_matrix(b) after loop for b */
#define update_matrix(b) do {payoff_matrix[1][0]=b;} while(0)

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
    for (i = 0; i < STRANUM; i++) stra_cnt[i] = 0;
    for (i = 0; i < SIZE; i++) stra_cnt[stra[i] = randi(2)]++;
}

double payoff(int x) {
    int i;
    double pay = 0;
    for (i = 0; i < NEINUM; i++)
        pay += payoff_matrix[stra[x]][stra[net[x][i]]];
    return pay;
}

void update_stra(int x) {
    int y = net[x][(int)randi(NEINUM)];
    /* Update strategy, if same then exit */
    if (stra[x] == stra[y]) return ;
    if (randf() < 1/(1+exp((payoff(x) - payoff(y))/K))) {
        stra_cnt[stra[x]]--;
        stra[x] = stra[y];
        stra_cnt[stra[x]]++;
    }
}

/* Snapshot function temporarily not used */
void snapshot(int step) {
    int x;
    sprintf(_file_path, "./%s/%s", OUTPUT_DIR, SNAPSHOT_DIR);
    make_dir(_file_path);
    sprintf(_file_path, "%s/%s/step=%d.%s", OUTPUT_DIR, SNAPSHOT_DIR, step+1, DATA_FN_EXT);
    FILE *fp = fopen(_file_path, "w");
    if (!fp) die(ENOENT, "Cannot create file for snapshot.");
    for (x = 0; x < SIZE; x++)
        fprintf(fp, "%c%s", STRA_STR[stra[x]], (x+1)%L?VALUE_SEP:"\n");
    fclose(fp);
}

void run() {
    /* Allocate net, stra */
    net = (_neigh *)malloc(SIZE * sizeof(_neigh));
    stra = (char *)malloc(SIZE * sizeof(char));
    if (!net || !stra) die(ENOMEM, "L is too large, no enough space.");
    FILE *fp_avg, *fp_steps;
    if (AVG_FN[0]) {
        sprintf(_file_path, "%s/%s_b=%g-%g.%s", OUTPUT_DIR, AVG_FN, bs, be, DATA_FN_EXT);
        fp_avg = fopen(_file_path, "w");
    }
    if (STEPS_DIR[0]) {
        sprintf(_file_path, "./%s/%s", OUTPUT_DIR, STEPS_DIR);
        make_dir(_file_path);
    }
    
    if (SEED != -1) sgenrand(SEED);
    else sgenrand((unsigned int)time(NULL));
    prod_neighbors();
    int x, step, rec, tac[TRY_TIMES], try_time;   /* total amount C in all REC steps */
    double fc, afc[TRY_TIMES];          /* Frequency of C */
    for (b = bs; b < be || fabs(b-be) < 1e-7; b += bi) {
        if (STEPS_DIR[0]) {
            sprintf(_file_path, "%s/%s/b=%g.%s", OUTPUT_DIR, STEPS_DIR, b, DATA_FN_EXT);
            fp_steps = fopen(_file_path, "w");
        }
        update_matrix(b);
        for (try_time = 0; try_time < TRY_TIMES; try_time ++) {
            init_randomly();
            tac[try_time] = rec = 0;
            for (step = 0; step < MC_STEPS; step++) {
                for (x = 0; x < SIZE; x ++)
                    update_stra((int)randi(SIZE));
                
                fc = stra_cnt[0]/(double)SIZE;
                if (step > MC_STEPS-AVG_STEPS-1) tac[try_time] += stra_cnt[0], rec++;
                
                if (STEPS_DIR[0]) {
                    if (TRY_TIMES > 1) fprintf(fp_steps, "%d%s", try_time, VALUE_SEP);
                	fprintf(fp_steps, "%d%s%g\n", step, VALUE_SEP, fc);
                }
            	
                if (REFRESH_FREQ && step % REFRESH_FREQ == 0) {
                    if (TRY_TIMES > 1) printf("\r [%d] Step: %d  C: %5.2lf%%          ", try_time, step, fc*100);
                    else printf("\r Step: %d  C: %5.2lf%%          ", step, fc*100);
                }
                
            	if ((!stra_cnt[0]) || (!stra_cnt[1])) {
            	    if (!rec) afc[try_time] = stra_cnt[0] ? 1 : 0;
                    break;
                };
            }
            if (rec) afc[try_time] = (double)tac[try_time]/(SIZE * rec);
            
            if (TRY_TIMES > 1) 
                printf("\r [%d] b: %g  avg_C: %5.2lf%%               \n", try_time, b, afc[try_time]*100);
            else printf("\r b: %g  avg_C: %5.2lf%%               \n", b, afc[try_time]*100);
        }
        if (STEPS_DIR[0]) fclose(fp_steps);
        double aafc=0;
        for (x = 0; x < TRY_TIMES; x ++)
            aafc += afc[x];
        aafc /= TRY_TIMES;
        if (AVG_FN[0]) fprintf(fp_avg, "%g%s%g\n", b, VALUE_SEP, aafc);
        if (TRY_TIMES > 1) printf("[%d tries]b: %g  avg_C: %5.2lf%%\n\n", TRY_TIMES, b, aafc*100);
    }
    if (AVG_FN[0]) fclose(fp_avg);
}

char *get_progname(char ** const argv) {
    char *slash_pos = strrchr(argv[0], '/');
    if (!slash_pos) slash_pos = strrchr(argv[0], '\\');
    return slash_pos ? slash_pos+1 : argv[0];
}

int main(int argc, char **argv) {
    progname = get_progname(argv);
    read_parameters(argc, argv);
    if (!STEPS_DIR[0] && !AVG_FN[0])
        fprintf(stderr, "[WARNING] Your output data will not be stored.\n");
    else save_metadata();
    run();
    return 0;
}

