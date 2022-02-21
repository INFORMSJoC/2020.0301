"""
Auth:   Wei WU.
Date:   2021/2/3
Desp:   Script for generating all results from paper.
        Run algorithms to solve binary integer programming (BIP) problems under the min-max regret criterion.
        For the following four special cases, the instances can be read directly:
            - the knapsack problem (KP),
            - the multidimensional knapsack problem (MKP),
            - the generalized assignment problem (GAP),
            - the set covering problem (SCP).
        Five types of algorithms are implemented:
            - fixed scenario algorithm (fix),
            - branch-and-cut algorithm (bc);
            - dual substitution algorithm (ds);
            - iterated dual substitution algorithm with best-scenario constraints (ids-b);
            - iterated dual substitution algorithm with Hamming-distance constraints (ids-h).
"""

from mmrbipy import Model
import os
import argparse

# path for results and data
DAT_PATH = '../data/'
RES_PATH = '../results/'

# define all the problem and algorithm names
ALGOS = ('fix','ds','bc','ids-b','ids-h',)
PROBS = ('KP','MKP','SCP','GAP',)

def run_all(prob: str='', algo: str='', timelimit: int=3600):
    """Run all test for all benchmark instances.

    Args:
        prob (str, optional): Problem name to test all. Defaults to '' (test all problems).
        algo (str, optional): Algorithm to run. Defaults to '' (test all algorithms).
        timelimit (int, optional): Time limit for each run. Defaults to 3600 (seconds).

    Raises:
        Exception: '../data/' direct not found.
    """
    if prob:
        prob = prob.upper()
    if algo:
        algo = algo.lower()
    probs = (prob,) if prob in PROBS else PROBS
    algos = (algo,) if algo in ALGOS else ALGOS

    if not os.path.exists(DAT_PATH):
        raise Exception("{} directory does not exist.".format(DAT_PATH))

    if not os.path.exists(RES_PATH):
        os.makedirs(RES_PATH)

    for prob in probs:
        path = DAT_PATH + prob
        inss = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        inss.sort()
        for algo in algos:
            wfile = "{}rst_{}_{}.csv".format(RES_PATH, prob, algo)
            with open (wfile, 'w') as f:
                f.write("instance,objective value,time to find best soltuion")
                if algo == 'bc':
                    f.write(",lower bound")
                if algo[:3] == 'ids':
                    f.write(",number of total iterations,number of iterations to find best soltuion")
                f.write("\n")
            for ins in inss:
                mod = Model(problem=prob, filename=os.path.join(path, ins))
                mod.solve(algorithm=algo, timelimit=timelimit)
                with open (wfile, 'a') as f:
                    f.write("{},{},{:.2f}".format(ins,mod.objval, mod.ttb))
                    if algo == 'bc':
                        f.write(",{}".format(mod.lb))
                    if algo[:3] == 'ids':
                        f.write(",{},{}".format(mod.num_ite, mod.bes_ite))
                    f.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solves all instances')
    # optional arguments
    parser.add_argument(
        '-p',
        '--problem',
        choices=['kp', 'mkp', 'gap', 'scp', 'bip'],
        help='Problem name among [kp, mkp, gap, scp, bip]')
    parser.add_argument(
        '-a',
        '--algorithm',
        choices=['fix', 'bc', 'ds', 'ids-h', 'ids-b'],
        help='Algorithm name among [fix, bc, ds, ids-h, ids-b]')

    parser.add_argument(
        '-t',
        '--timelimit',
        type=int,
        default=3600,
        help='Time limit for solving the each instance')

    args = parser.parse_args()
    run_all(prob=args.problem, algo=args.algorithm, timelimit=args.timelimit)
    exit(0)