#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Random function */
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
/* End of random */














#define RANDOMIZE   3145215
#define L           200
#define K           0.1
#define NEINUM      4
#define SIZE        (L * L)
#define MC_STEPS    10000 /*50000 steps for general*/
#define REC_STEPS   5000
#define REFRESH_FRE 100     /* Control the frequency of refresh screen, 1 for all the time */
#define TRY_TIME    10
#define b           1.00
#define OUTER       0.16
#define IN          0.3
/* green area:OUTER=0.16 and IN=0.3*/

int net[SIZE][NEINUM];   /* Players neighbor relations */
int stra[SIZE], pay[SIZE];

/* Properties counter */
int stra_cnt[2] = { 0 };    /*strategy matrix: 0 for C and 1 for D with its initialized value: 0*/
int pay_cnt[2] = { 0 };


/* Payoff matrix and its update */
double payoff_matrix[2][2] = {
    {1, 0},
    {0, 0}
};
/* Call update_matrix(b) after loop for b */
#define update_matrix(b) payoff_matrix[1][0]=b
/*payoff_matrix[1][0]表示d和c交手时,存在b的收益*/


/* Construct Clockwise direction */ 
char CW[] = {1, 3, 0, 2};
char CC[] = {2, 0, 3, 1};
char OP[] = {3, 2, 1, 0};

void prod_neighbors() {
    int i, j, x;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++)
        {
            x = i*L + j; //convert lattice into one demension
            net[x][0] = i*L + ((j - 1 + L) % L);  /* left */
            net[x][1] = ((i - 1 + L) % L)*L + j;  /* up */
            net[x][2] = ((i + 1) % L)*L + j;      /* down */
            net[x][3] = i*L + ((j + 1) % L);      /* right */
            /* its four neighbor */
        }
}

// Construct Strategies Matrix
#define ST_MAT_ROW  2
#define ST_MAT_COL  2
const char ST_MAT[ST_MAT_ROW][ST_MAT_COL] = {
    {1, 0},
    {0, 1}
};

/* === Payoff matrix distribution ===
         8  11
         |   |
     6-- 7--10--12            op
     |   |   |   |
 5-- 4-- 0-- 9--13--14     cw    cc
     |   |   |   |
     3-- 1--16--15            di
         |   |
         2  17
      x = 0, y = 9
      di is the direction below x,
      and the direction above y.
Payoffs:

     7  10              4   5
     |   |              |   |
 4-- 0-- 9--13 ==>  3-- 0-- 1-- 6
     |   |              |   |
     1  16              2   7
*/

/* Neighbor map */
const int NM[8][4] = {  
    { 4,  7,  1,  9},
    { 0, 10, 16, 13},
    { 3,  0,  2, 16},
    { 5,  6,  3,  0},
    { 6,  8,  0, 10},
    { 7, 11,  9, 12},
    { 9, 12, 15, 14},
    { 1,  9, 17, 15},
};
/*  Construct position array, input array, 
    current player id and init direction */
void con_pos(int *pos, int x, char di) {
    char cw = CW[di], op = OP[di];
    /* di is original direction. */
    pos[0] = x;
    pos[1] = net[x][di];
    pos[2] = net[pos[1]][di];
    pos[3] = net[pos[1]][cw];
    pos[4] = net[x][cw];
    pos[5] = net[pos[4]][cw];
    pos[6] = net[pos[4]][op];
    pos[7] = net[x][op];
    pos[8] = net[pos[7]][op];
}

/* payoff caculation */
void cal_payoffs(double *pap, int *pa) {
    int i, j;
    for (i = 0; i < 18; i++) {
        pap[i] = 0;
        for (j = 0; j < NEINUM; j++)
            pap[i] += payoff_matrix[
                stra[pa[i]]][stra[net[pa[i]][j]]];
    }
}

/*test:print the payoff matrix*/
void cal_avgs(double *avgs, double *pap) {
    int i, j;
    for (i = 0; i < 8; i++) {
        avgs[i] = 0;
        for (j = 0; j < NEINUM; j++)
            avgs[i] += pap[NM[i][j]];
    }
}

/* new_str due to community(the same strategy) */
double new_payoff(int *pa, double *pap, double *avgs, int pos) {
    int i, x = pa[pos];;
    double new_pay=pap[pos];
    /* as to itself */
    if(new_pay > avgs[pos]) {
        /*change from c to a rc or from p to a rp */
        for (i = 0; i < NEINUM; i++)
            if (stra[x] == stra[net[x][i]])
                new_pay -= OUTER;
        if(pay[x] == 0) {
            pay_cnt[pay[x]]--;
            pay[x] = 1;
            pay_cnt[pay[x]]++;
        }
    }
    else {
        if(pay[x] == 1) {
            pay_cnt[pay[x]]--;
            pay[x] = 0;
            pay_cnt[pay[x]]++;
        }
    }
    for (i = 0; i < NEINUM; i++)
        if(stra[net[x][i]] == stra[x])
            if(pap[NM[pos>0][i]]>avgs[NM[pos>0][i]])
                new_pay += IN;
    return new_pay;
}

/* strategy chosen */
void update_stra(int x) {
    int i=randi(4);
    int y=net[x][i];
    if (stra[x] == stra[y]) return;
    int pa[18];      /* Construct position array */
    double pap[18], avgs[8];
    con_pos(pa, x, CW[i]);
    con_pos(pa+9, y, CC[i]);
    cal_payoffs(pap, pa);
    cal_avgs(avgs, pap);
    double x_pay = new_payoff(pa, pap, avgs, 0), y_pay = new_payoff(pa, pap, avgs, 9);
    /* Update strategy, if same then exit */
    if (randf() < 1 / (1 + exp((x_pay - y_pay) / K))) {
        stra_cnt[stra[x]]--;
        stra[x] = stra[y]; //learn the basic strategy of y
        stra_cnt[stra[x]]++;
    }
}


void init_by_st_mat() {
    stra_cnt[0] = stra_cnt[1] = pay_cnt[0] = pay_cnt[1] = 0;
    int i, j, block_row = L/ST_MAT_ROW, block_col = L/ST_MAT_COL;
    for (i = 0; i < L; i++)
        for (j = 0; j < L; j++) {
            stra_cnt[stra[i*L+j] = ST_MAT[i/block_row][j/block_col]]++;
            pay_cnt[pay[i*L+j] = 0]++;
        }
}

void init() {
    stra_cnt[0] = stra_cnt[1] = pay_cnt[0] = pay_cnt[1] = 0;
    int i, j; //i is the row;j is the line
    for (i = 0; i < SIZE; i++) {
        stra_cnt[stra[i] = randi(2)]++;  //first: there is no rc and rp
        pay_cnt[pay[i] = 0]++;
    }
}

void snapshot(int step) {
    char fpath[64];
    int x;
    sprintf(fpath, "snapshots/lattice_step=%d.txt", step+1);
    FILE *fp = fopen(fpath, "w");
    if (!fp) {
        fprintf(stderr, "Maybe folder (snapshots) is not created.\n");
        exit(-1);
    }
    for (x = 0; x < SIZE; x++)
        fprintf(fp, "%c%c%c", "OR"[pay[x]], "CD"[stra[x]], (x+1)%L?',':'\n');
}

int main() {
    sgenrand(RANDOMIZE);
    prod_neighbors();
    double fc;
    int x, step;
    update_matrix(b);//博弈矩阵
    init();
    /*
    for (x = 0; x < SIZE; x++)
        fprintf(fp_lattice_begin, "%d\n", stra[x]);
    */
    snapshot(-1);       // Before MCS
    int total=stra_cnt[0]+stra_cnt[1]+stra_cnt[2]+stra_cnt[3];
    for (step = 0; step < MC_STEPS; step++) {
        fc = (stra_cnt[0] +stra_cnt[2]) / (double)SIZE;
        if (step % REFRESH_FRE == 0)
            printf("\rStep: %d\t C: %lf%%\t  PAY: %lf%%         ", 
                step, stra_cnt[0]/(double)SIZE, pay_cnt[1]/(double)SIZE);
        
        for (x = 0; x < SIZE; x++)
        {
            update_stra((int)randi(SIZE));
            //printf("!");
            //return 0;
        }
        if (step == 0 || step == 5 || step == 10 || step == 50 || 
                step == 100)
            snapshot(step);
    }
    snapshot(step);
    putchar(10);
    return 0;
}

