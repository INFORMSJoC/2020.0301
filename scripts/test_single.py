"""
Auth:   Wei WU.
Date:   2021/2/3
Desp:   Script for generating results from paper. Single run for one instance unisng one algorithm.
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

import argparse
from mmrbipy import Model

def main():
    parser = argparse.ArgumentParser(description='Solves a MMR-BIP')
    # required arguments
    parser.add_argument(
        '-p',
        '--problem',
        choices=['kp', 'mkp', 'gap', 'scp', 'bip'],
        required=True,
        help='Problem name among [kp, mkp, gap, scp, bip]')
    parser.add_argument(
        '-i',
        '--instance',
        required=True,
        help='Filename for the problem instance')

    # optional arguments
    parser.add_argument(
        '-a',
        '--algorithm',
        choices=['fix', 'bc', 'ds', 'ids-h', 'ids-b'],
        default='ids-b',
        help='Algorithm name among [fix, bc, ds, ids-h, ids-b]')
    parser.add_argument(
        '-t',
        '--timelimit',
        type=int,
        default=3600,
        help='Time limit for solving the problem')
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Name of file to which to write the result")

    args = parser.parse_args()
    mod = Model(problem=args.problem, filename=args.instance)
    mod.solve(algorithm=args.algorithm, timelimit=args.timelimit)
    if args.output:
        mod.write(args.output)

    print("objective value: {}".format(mod.objval))
    print("time to best: {:.2f}".format(mod.ttb))
    if args.algorithm == 'bc':
        print("lower bound: {}".format(mod.lb))
    if args.algorithm[:3] == 'ids':
        print("#iteration: {}".format(mod.num_ite))
        print("best iteration: {}".format(mod.bes_ite))
    
if __name__ == "__main__":
    main()
    exit(0)