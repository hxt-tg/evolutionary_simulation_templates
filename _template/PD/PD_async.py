import networkx as nx
import numpy as np
from collections import namedtuple
import os

Range = namedtuple("Range", ["start", "end", "interval"])

# Global settings
REC_SNAPSHOT  = 1             # Not equal to 0 if record snapshots
REC_STEPS     = 1             # Not equal to 0 if record data per step
RANDSEED      = 3145215       # Fixed random seed, use time if -1
# Network settings
NET_TYPE      = "GRID"        # Network type
L             = 100           # Side length of a network
SIZE          = (L * L)       # Total number of nodes
N_NEIGH       = 4             # Number of neighbors for every node
# Monte Carlo settings
K             = 0.1           # Fermi coefficient
MC_STEPS      = 10000         # Number of MC steps
AVG_STEPS     = 5000          # Number of steps at the end for averaging
REFRESH_FREQ  = 100           # Control the frequency of refresh screen, 1 for all the time
TRY_TIMES     = 1             # Average result in tried times

# Strategies
_C = 0
_D = 1
N_STRA = 2

class Game:
    def __init__(self, b_range: Range):
        self.b_range = b_range
        self.b = self.b_range.start
        self.stra_cnt = [0] * N_STRA
        self._rand = -1 if RANDSEED == -1 else RANDSEED
        self._payoff_matrix = np.array([
            [1, 0],
            [0, 0]
        ])
        if self._rand != -1:
            np.random.seed(self._rand)

    def start(self):
        try:
            os.makedirs("output")
            if REC_STEPS: os.makedirs("output/steps")
            if REC_SNAPSHOT: os.makedirs("output/snapshots")
        except FileExistsError:
            pass

        print("Construct network ... ", end='')
        if NET_TYPE == 'GRID':
            self.net = nx.grid_2d_graph(L, L, periodic=True)
            node_map = {n: n[0]*L+n[1] for n in self.net}
            self.net = nx.relabel_nodes(self.net, node_map)
        else:
            raise ValueError("Unspecified net structure.")
        for i in range(SIZE):
            self.net.nodes[i]["stra"] = _C
            self.stra_cnt[self.net.nodes[i]["stra"]] += 1
        print("Done")
        self._save_meta()

        fp_aver = open("output/average.csv", "w")
        fp_aver.write("b,f_C,f_D\n")
        for self.b in np.arange(*self.b_range):
            self._payoff_matrix[1, 0] = self.b
            self._run(fp_aver)
        fp_aver.close()

    def _run(self, fp_aver):
        if REC_STEPS:
            fp_steps = open(f"output/steps/b={self.b:.2f}.csv", "w")
            if TRY_TIMES > 1: fp_steps.write("try,")
            fp_steps.write("step,f_C,f_D\n")
        f = [0]*N_STRA
        af = [[0]*N_STRA for _ in range(TRY_TIMES)]
        for try_time in range(TRY_TIMES):
            ta = [0] * N_STRA
            rec = 0
            self._init_randomly()
            self.snapshot(try_time, -1)
            for step in range(MC_STEPS):
                for x in range(SIZE):
                    self._update_stra(np.random.randint(SIZE))
                for i in range(N_STRA):
                    f[i] = self.stra_cnt[i]/SIZE

                if step > MC_STEPS-REC_STEPS - 1:
                    for i in range(N_STRA):
                        ta[i] += self.stra_cnt[i]
                    rec += 1

                if TRY_TIMES > 1: fp_steps.write(f"{try_time},")
                fp_steps.write(f"{step},{f[_C]},{f[_D]}\n")

                if REFRESH_FREQ and step % REFRESH_FREQ == 0 :
                    print("\r" if TRY_TIMES == 1 else f"\r [{try_time}]", end='')
                    print(f" Step: {step:6d}  f: [{f[_C]*100:6.2f}%,{f[_D]*100:6.2f}%]", end='')
                if self.stra_cnt[_C] == SIZE or self.stra_cnt[_D] == SIZE:
                    if step < MC_STEPS-AVG_STEPS:
                        for i in range(N_STRA):
                            af[try_time][i] = 1 if self.stra_cnt[i] else 0
                    step += 1
                    break
            if rec:
                for i in range(N_STRA):
                    af[try_time][i] = ta[i]/(SIZE * rec)

            print("\r" if TRY_TIMES == 1 else f"\r [{try_time}]", end='')
            print(f" b: {self.b:.2f}   avg: [{af[try_time][_C]*100:6.2f}%,{af[try_time][_D]*100:6.2f}%]       ")

        if REC_STEPS: fp_steps.close()
        aaf = [0] * N_STRA
        for i in range(N_STRA):
            for x in range(TRY_TIMES):
                aaf[i] += af[x][i]
            aaf[i] /= TRY_TIMES
        fp_aver.write(f"{self.b:.2f},{aaf[_C]:f},{aaf[_D]:f}\n")
        if TRY_TIMES > 1: print(f"[{TRY_TIMES} tries] b: {self.b:g}  avg: [{aaf[_C]*100:6.2f}%,{aaf[_D]*100:6.2f}%]\n")

    def _payoff(self, x):
        # Isolated node
        if self.net.degree(x) == 0: return 0
        pay = 0
        for n in self.net.neighbors(x):
            pay += self._payoff_matrix[self.net.nodes[x]['stra']][self.net.nodes[n]['stra']]
        return pay

    def _update_stra(self, x):
        # Isolated node
        if self.net.degree(x) == 0: return 0
        y = np.random.choice(list(self.net.neighbors(x)))
        # Update strategy, if same then exit
        if self.net.nodes[x]['stra'] == self.net.nodes[y]['stra']: return
        if np.random.rand() < 1/(1+np.exp((self._payoff(x) - self._payoff(y))/K)):
            self.stra_cnt[self.net.nodes[x]['stra']] -= 1
            self.net.nodes[x]['stra'] = self.net.nodes[y]['stra']
            self.stra_cnt[self.net.nodes[y]['stra']] += 1

    def snapshot(self, try_time, step):
        if not REC_SNAPSHOT: return
        if TRY_TIMES > 1: filepath = f"output/snapshots/b={self.b:.2f}_step={step+1}[{try_time}].csv"
        else: filepath = f"output/snapshots/b={self.b:.2f}_step={step+1}.csv"
        with open(filepath, "w") as fp:
            for x in range(SIZE):
                fp.write(f"{self.net.nodes[x]['stra']:d}"+(','if (x+1)%L else ''))

    def _init_randomly(self):
        self.stra_cnt = [0] * N_STRA
        for i in range(SIZE):
            self.net.nodes[i]['stra'] = np.random.randint(N_STRA)
            self.stra_cnt[self.net.nodes[i]['stra']] += 1

    def _save_meta(self):
        with open("output/meta.csv", "w") as f:
            f.write(f"L,{L},Side length of square lattice\n")
            f.write(f"SIZE,{SIZE},Number of agents on network\n")
            f.write(f"N_NEIGH,{N_NEIGH},Number of neighbors for every agent\n")
            f.write(f"K,{K},Fermi coefficient\n")
            f.write(f"_rand,{self._rand},Random seed\n")
            f.write(f"MC_STEPS,{MC_STEPS},Monte Carlo steps\n")
            f.write(f"AVG_STEPS,{AVG_STEPS},Recorded MC steps at the end\n")
            f.write(f"TRY_TIMES,{TRY_TIMES},Try times per parameter\n")
            f.write(f"b={self.b_range.start:g}~"
                    f"{self.b_range.end:g} with interval "
                    f"{self.b_range.interval:g}.\n")

g = Game(Range(1, 1.5, 0.1))
g.start()
