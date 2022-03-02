[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An iterated dual substitution approach for min-max binary integer programming

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported in the paper _An iterated dual substitution approach for min-max binary integer programming_ by W. Wu, M. Iori, S. Martello, and M. Yagiura. 

**Important: This code is being developed on an on-going basis at 
https://github.com/ebeleta/iDS. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite this software, please cite the software itself, using the following DOI.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6227712.svg)](https://doi.org/10.5281/zenodo.6227712)

Below is the BibTex for citing this version of the code.

```
@article{iDS2022,
  author =        {W. Wu, M. Iori, S. Martello, and M. Yagiura},
  publisher =     {INFORMS Journal on Computing},
  title =         {{mmrbipy} Version v2020.0301},
  year =          {2022},
  doi =           {10.5281/zenodo.6227712},
  url =           {https://github.com/INFORMSJoC/2020.0301},
}  
```

## Installation

In a virtual environment with Python 3.6+, mmrbipy can be installed via

```bash
pip install mmrbipy
```

## Using mmrbipy

With a compatible instance file, mmrbipy solves the MMR-BIP from a Python script:

```python
from mmrbipy import Model

# Generate a model from an instance file
mod = Model(problem='kp', filename='../data/KP/1-70-01-45-20')

# Solve the model by using iDS algorithm with best-scenario constraints
mod.solve(algorithm='ids-b', timelimit=100)

# Print results
print("objective value: {}".format(mod.objval))
print("time to best: {:.2f}".format(mod.ttb))

# Write results to a text file
mod.write("result.txt")
```
## Model
To solve the MMR-BIP, mmrbipy provides four types of instance formats:

- min-max regret knapsack problem (*kp*)
- min-max regret multidimensional knapsack problem (*mkp*)
- min-max regret set covering problem (*scp*)
- min-max regret generalized assignment problem (*gap*)

See [data](data) directory for the details of each type.

### Set problem type in constructor of _Model_ class
```python
# Generate a model from an instance file
mod = Model(problem='kp', filename='../data/KP/1-70-01-45-20')
```

_Note: Benchmark instances for_

- _min-max regret knapsack problem_,
- _min-max regret multidimensional knapsack problem_,
- _min-max regret set covering problem_,
- _min-max regret generalized assignment problem_

_are available in [data](data) directory._

## Algorithms

To solve the MMR-BIP, mmrbipy provides five algorithms:
- fixed scenario algorithm (*fix*),
- branch-and-cut algorithm (*bc*),
- dual substitution algorithm (*ds*),
- iterated dual substitution algorithm with best-scenario constraints (*ids-b*),
- iterated dual substitution algorithm with Hamming-distance constraints (*ids-h*).

### Set algorithm type in _solve_ function
```python
# Solve the model by using iDS algorithm with best-scenario constraints
mod.solve(algorithm='ids-b', timelimit=100)
```

_Note: For this project, we use [gurobipy](https://pypi.org/project/gurobipy/) as a solver to solve mixed integer programming problems._

## Results

Detailed results for all the tested instances can be found in [result](results/results.pdf) pdf.
We used Gurobi Optimizer version 8.1 to solve mixed integer programming problems.

## Ongoing development

This code is being developed on an on-going basis at the author's
[Github site](https://github.com/ebeleta/iDS).

## Support

For support in using this software, submit an
[issue](https://github.com/ebeleta/iDS/issues/new).
