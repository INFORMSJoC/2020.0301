# Script for generating results from paper

Run `test_all.py` or `test_single.py` to replicate the results for an instance of the following problem  
- min-max binary integer programming problem (*bip*)
- min-max regret knapsack problem (*kp*)
- min-max regret multidimensional knapsack problem (*mkp*)
- min-max regret set covering problem (*scp*)
- min-max regret generalized assignment problem (*gap*)

by algorithms 
- fixed scenario algorithm (*fix*);
- branch-and-cut algorithm (*bc*);
- dual substitution algorithm (*ds*);
- iterated dual substitution algorithm with best-scenario constraints (*ids-b*);
- iterated dual substitution algorithm with Hamming-distance constraints (*ids-h*).

## Examples
### Test single
Run iterated dual substitution algorithm with best-scenario constraints to solve an MMR-KP instance with a 10-second time limit:

```bash
python test_single.py -i ../data/KP/1-70-01-45-20 -p kp -a ids-b -t 10 -o result.sol
```
The results will be written to `result.sol`.
### Test all
Run branch-and-cut algorithm to solve all the instances of the MMR-SCP with a 600-second time limit:

```bash
python test_all.py -p scp -a bc -t 600
```
The results will be written to `../results/rst_SCP_bc.csv`.
___

Run all the algorithms to solve all instances of all the problems with a 3600-second time limit:

```bash
python test_all.py
```
All the results will be written to `../results/` as csv files.

## Requirements
- _argparse_ and _mmrbipy_ modules are installed.
- For running `test_all.py`, `../data/` directory should exist.