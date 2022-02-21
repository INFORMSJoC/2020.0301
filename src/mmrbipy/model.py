"""
Module for solving min-max regret binary programming problem from Wu et al. (2022).
Copyright (c) 2022 W. Wu, M. Iori, S. Martello & M Yagiura 
Released under the MIT license
https://opensource.org/licenses/mit-license.php
"""

from gurobipy import GRB
import gurobipy as gp
import time

# Constant for rounding
EPSN = 0.0001

def mipsol(model, where):
    """Callback to obtain runtime when the incumbent is updated"""
    if where == GRB.Callback.MIPSOL:
        model._ttb = model.cbGet(GRB.Callback.RUNTIME)

def lazy_cut(model, where):
    """Callback to add a cut for branch-and-cut framework"""
    # Excute the function when an incumbent is found
    if where != GRB.Callback.MIPSOL:
        return
    # Obtain the incumbent solution
    x_sol = model.cbGetSolution(model._x)
    r_sol = model.cbGetSolution(model._r)
    # Obtain the Model object
    obj = model._obj
    
    # Generate the worst-case scenario for the solution
    N = obj.c_low.keys()
    if obj.sense == GRB.MAXIMIZE:
        c_wst = {j: (obj.c_low[j] - obj.c_upp[j]) * x_sol[j] + obj.c_upp[j] for j in N}
    else:
        c_wst = {j: (obj.c_upp[j] - obj.c_low[j]) * x_sol[j] + obj.c_low[j] for j in N}

    # Solve BIP under the worst-case scenario
    remain_time = obj.timelimit-model.cbGet(GRB.Callback.RUNTIME)
    res = obj.solve_bip(c=c_wst, timelimit=remain_time)

    if obj.sense == GRB.MAXIMIZE and res['obj_bd'] > r_sol + EPSN:
        # Add a cut
        rhs = gp.LinExpr()
        for j in N:
            if res['sol'][j] > 0.5:
                rhs += obj.c_upp[j]
                if obj.c_upp[j] - obj.c_low[j] > EPSN:
                    rhs -= (obj.c_upp[j] - obj.c_low[j]) * model._x[j]
        model.cbLazy(model._r >= rhs)
    elif obj.sense == GRB.MINIMIZE and res['obj_bd'] < r_sol - EPSN:
        # Add a cut
        rhs = gp.LinExpr()
        for j in N:
            if res['sol'][j] > 0.5:
                rhs += obj.c_low[j]
                if obj.c_upp[j] - obj.c_low[j] > EPSN:
                    rhs += (obj.c_upp[j] - obj.c_low[j]) * model._x[j]
        model.cbLazy(model._r <= rhs)
    else:
        # If the incumbent is feasible, stamp the runtime
        model._ttb = model.cbGet(GRB.Callback.RUNTIME)

class Model():
    """A class to model and solve a BIP problem instance.
    Attributes:
        filename (str): instance file name
        problem (str): problem name

        sol (dict): incumbent solution
        objval (int): objective value
        lb (int): lower bound of the objective value
        ttb (float): time to best (runtime to obtain the best solution)
        bes_ite (int): iteration to best (for iDS)
        num_ite (int): total iteration numbers (for iDS)
    """
    
    def __init__(self, problem: str, filename: str):
        """[summary]

        Args:
            problem (str): problem type 
                           (kp: knapsack problem, mkp: multidimensional knapsack problem
                            gap: generalized assignment problem, scp: set covering problem
                            bip: binary programming problem)
            filename (str): instance filename
        """
        # Read the data file
        with open(filename) as file:
            lines = file.readlines()
        lines = [line.replace("\n", "").replace("\r", "").rstrip().replace("\t", " ") for line in lines]
        
        # Validate problem type
        problem = problem.lower()
        if problem not in {'kp', 'mkp', 'gap', 'scp', 'bip'}:
            raise Exception("Problem type is not correct")
        # Store the instance date in class object
        self.a, self.b, self.c_upp, self.c_low = dict(), dict(), dict(), dict()
        # Read instance and translate it to a standard BIP
        if problem == 'kp':
            # Read KP benchmark instance
            self.sense = GRB.MAXIMIZE
            _, self.b[0] = int(lines.pop(0)), int(lines.pop(0))
            for j,v in enumerate(lines.pop(0).split()):
                self.a[0,j] = int(v)
            for j,v in enumerate(lines.pop(0).split()):
                self.c_low[j] = int(v)
            for j,v in enumerate(lines.pop(0).split()):
                self.c_upp[j] = int(v)
            self.m = 1
        elif problem == 'mkp':
            # Read MKP benchmark instance
            self.sense = GRB.MAXIMIZE
            self.m, _ = [int(item) for item in lines.pop(0).split()]
            for j,v in enumerate(lines.pop(0).split()):
                self.c_low[j] = int(v)
            for j,v in enumerate(lines.pop(0).split()):
                self.c_upp[j] = int(v)
            for i in range(self.m):
                for j,v in enumerate(lines.pop(0).split()):
                    self.a[i,j] = int(v)
            for i,v in enumerate(lines.pop(0).split()):
                self.b[i] = int(v)
        elif problem == 'gap':
            # Read GAP benchmark instance
            self.sense = GRB.MINIMIZE
            m, n = int(lines.pop(0)), int(lines.pop(0))
            for i in range(m):
                for j,v in enumerate(lines.pop(0).split()):
                    self.c_low[i,j] = int(v)
            for i in range(m):
                for j,v in enumerate(lines.pop(0).split()):
                    self.c_upp[i,j] = int(v)
            for i in range(m+2*n):
                for (ii, jj) in self.c_upp:
                    self.a[i,(ii, jj)] = 0
            for i in range(m):
                for j,v in enumerate(lines.pop(0).split()):
                    self.a[i,(i,j)] = int(v)
            for i,v in enumerate(lines.pop(0).split()):
                self.b[i] = int(v)
            for j in range(n):
                for i in range(m):
                    self.a[m+j,(i, j)] = 1
                    self.a[m+n+j,(i, j)] = -1
                self.b[m+j] = 1
                self.b[m+n+j] = -1
            self.m = m + 2 * n
        elif problem == 'scp':
            # Read SCP benchmark instance
            self.sense = GRB.MINIMIZE
            self.m, n = [int(item) for item in lines.pop(0).split()]
            for j in range(n):
                self.c_low[j], self.c_upp[j] = [int(item) for item in lines.pop(0).split()]
            self.a = {(i,j): 0 for i in range(self.m) for j in range(n)}
            for i in range(self.m):
                idx = [int(item) for item in lines.pop(0).split()]
                for k in range(idx[0]):
                    j = idx[k+1]
                    self.a[i,j] = -1
            self.b = {j: -1 for j in range (self.m)}
        elif problem == 'bip':
            # Read BIP instance
            sense = lines.pop(0)[:3].upper()
            if sense == 'MIN':
                self.sense = GRB.MINIMIZE
            elif sense == 'MAX':
                self.sense = GRB.MAXIMIZE
            self.m, _ = [int(item) for item in lines.pop(0).split()]
            for j,v in enumerate(lines.pop(0).split()):
                self.c_low[j] = int(v)
            for j,v in enumerate(lines.pop(0).split()):
                self.c_upp[j] = int(v)
            for i in range(self.m):
                for j,v in enumerate(lines.pop(0).split()):
                    self.a[i,j] = int(v)
            for i,v in enumerate(lines.pop(0).split()):
                self.b[i] = int(v)
    
    def solve(self, algorithm: str = 'fix', timelimit: int = 30):
        """Solve the problem

        Args:
            algorithm (str): algorithm name. Defaults fix.
            timelimit (int, optional): time limit. Defaults 30.
        """
        if algorithm not in {'fix', 'bc', 'ds', 'ids-h', 'ids-b'}:
            raise Exception("Algorithm type is not correct")
        self.algorithm = algorithm
        # Set time limit
        self.timelimit = timelimit
        # Set start time
        self.starttime = time.time()
        # Run algorithm 
        if algorithm == 'fix':
            self.algo_fix()
        elif algorithm == 'ds':
            self.algo_ds()
        elif algorithm == 'bc':
            self.algo_bc()
        elif algorithm[:3] == 'ids':
            self.algo_ids(cut_type=algorithm[-1])
        self.algo = algorithm

    def algo_ids(self, cut_type: str):
        """Solve the problem by iDS algorithms

        Args:
            cut_type (str): cut type using in iDS
        """
        # Initialize local variables
        time_remain, ite_num = self.timelimit, 0
        best_obj, best_ite, best_time, best_sol = sum(self.c_upp.values()), -1, None, None
        # Initialize the model as a DS model
        model, x = self.set_ds_model()

        while time_remain > 0:
            # Solve the current DS model
            model.Params.timeLimit = time_remain
            model.optimize(mipsol)
            # Log the iteration number
            ite_num += 1
            # If the model is infeasible, an exact solution is obtained.
            if model.status == GRB.INFEASIBLE: break
            # If no feasible solution found, terminate the approach
            if model.solCount <= 0: break
            # Compute the regret of the current solution
            sol = {j: 1 if x[j].x > 0.5 else 0 for j in x}
            regret = self.evaluate(sol)
            # Update the best if a smaller regret is obtained
            if regret < best_obj - EPSN:
                best_obj, best_sol, best_ite = regret, sol, ite_num
                best_time = self.timelimit - time_remain + model._ttb
                print("#ite = {},\tobj = {}\ttime = {:.3f}".format(best_ite, best_obj, best_time))
            # Add a constraint based on the current solution
            lhs = gp.LinExpr()
            if cut_type == 'h':
                # Cut using Hamming distance 1
                rhs = 1
                for j in x:
                    if x[j].x < 0.5:
                        lhs += x[j]
                    else:
                        lhs -= x[j]
                        rhs -= 1
                model.addConstr((lhs >= rhs), name='ITE_CON_'+str(ite_num))
            elif cut_type == 'b':
                # Cut using Best-Scenario Lemma
                rhs = 1
                for j in x:
                    if self.sense == GRB.MAXIMIZE:
                        if x[j].x < 0.5:
                            lhs += self.c_upp[j] * x[j]
                        else:
                            lhs += self.c_low[j] * x[j]
                            rhs += self.c_low[j]
                    else:
                        if x[j].x < 0.5:
                            lhs -= self.c_low[j] * x[j]
                        else:
                            lhs -= self.c_upp[j] * x[j]
                            rhs -= self.c_upp[j]
                model.addConstr((lhs >= rhs), name='ITE_CON_'+str(ite_num))
            # Update remaining time
            time_remain = self.timelimit - (time.time() - self.starttime)
        # Store results
        self.objval = int(best_obj + EPSN)
        self.ttb = best_time
        self.sol = best_sol
        self.bes_ite = best_ite
        self.num_ite = ite_num

    def algo_fix(self):
        """Solve BIP under median scenario and evaluate the obtained solution under its worst-case scenario
        """
        # Generate median scenario
        c = {j: float(self.c_upp[j] + self.c_low[j])/2 for j in self.c_upp}
        # Solve BIP under median scenario
        res = self.solve_bip(c=c, timelimit=self.timelimit)
        # Evaluate the obtained solution under its worst-case scenario
        if res:
            self.sol = res['sol']
            self.ttb = res['ttb']
            self.objval = self.evaluate(res['sol'])

    def algo_ds(self):
        """Solve MMR-BIP by dual substitution approach
        """
        # Generate GUROBI model
        model, x = self.set_ds_model(timelimit=self.timelimit)
        # Use GUROBI to optimize
        model.optimize(mipsol)
        # Store results
        if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
            self.sol = {j: 1 if x[j].x > 0.5 else 0 for j in x}
            self.ttb = model._ttb
            # Evaluate the obtained solution under its worst-case scenario
            self.objval = self.evaluate(self.sol)

    def algo_bc(self):
        """Solve BIP by branch-and-cut approach
        """
        N, M = self.c_low.keys(), range(self.m)
        # Build model
        model = gp.Model('bc')
        # Define variables
        x = {j: model.addVar(vtype=GRB.BINARY, name='x[{}]'.format(j)) for j in N}
        r = model.addVar(vtype=GRB.CONTINUOUS, name='lambda', ub=sum(self.c_upp.values()))
        model.update()
        # Objective function
        if self.sense == GRB.MAXIMIZE:
            model.setObjective(r - gp.quicksum(self.c_low[j] * x[j] for j in N), GRB.MINIMIZE)
        else:
            model.setObjective(gp.quicksum(self.c_upp[j] * x[j] for j in N) - r, GRB.MINIMIZE)
        # Add constraints
        model.addConstrs((gp.quicksum(self.a[i,j] * x[j] for j in N) <= self.b[i] for i in M), name='CAP_CON')
        # Store objects for callback
        model._x = x
        model._r = r
        model._obj = self
        # Parameter setting
        model.Params.outputFlag = False
        model.Params.threads = 1
        model.Params.lazyConstraints = 1
        model.Params.timeLimit = self.timelimit
        model.Params.MIPGap = 0.0
        # Use GUROBI to optimize
        model.optimize(lazy_cut)
        # Store results
        if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
            self.lb = int(model.objBound + 1 - EPSN)
            self.objval = int(model.objVal + EPSN)
            self.ttb = model._ttb
            self.sol = {j: 1 if x[j].x > 0.5 else 0 for j in x}

    def evaluate(self, sol: dict) -> int:
        """Calculate regret for a solution

        Args:
            sol (dict): given solution

        Returns:
            int: regret value
        """
        # Prepare worst-case scenario
        if self.sense == GRB.MAXIMIZE:
            c_wst = {j: self.c_low[j] if sol[j] > 0.5 else self.c_upp[j] for j in self.c_upp}
        else:
            c_wst = {j: self.c_low[j] if sol[j] < 0.5 else self.c_upp[j] for j in self.c_upp}
        # Calculate the regret for the given solution
        val_sol = sum(c_wst[j] * sol[j] for j in self.c_upp)
        res_bst = self.solve_bip(c=c_wst)
        if self.sense == GRB.MAXIMIZE:
            regret = int(res_bst['obj_bd'] - val_sol + 0.5)
        else:
            regret = int(val_sol - res_bst['obj_bd'] + 0.5)
        return regret
    
    def solve_bip(self, c: dict, timelimit: int=0) -> dict:
        """Solve BIP by using GUROBI

        Args:
            c (dict): coefficients of objective function
            timelimit (int, optional): time limit. No limits if leave it as default value 0.

        Returns:
            dict: results including
                        solution (key:'sol', [dict]),
                        objective value (key:'obj_val', [float]),
                        bound of objective value (key:'obj_bd', [float]),
                        time to best (key:'ttb', [float]).
        """
        # Build model
        model = gp.Model('bip')
        # Define variables
        x = {j: model.addVar(vtype=GRB.BINARY, name='x[{}]'.format(j)) for j in c}
        model.update()
        # Objective function
        model.setObjective(gp.quicksum(c[j] * x[j] for j in x), self.sense)
        # Add constraints
        model.addConstrs((gp.quicksum(self.a[i,j] * x[j] for j in x) <= self.b[i] for i in range(self.m)), name='CON')
        # Parameter setting
        model.Params.outputFlag = False
        model.Params.threads = 1
        model.Params.MIPGap = 0.0
        model.Params.lazyConstraints = 1
        if timelimit:
            model.Params.timeLimit = timelimit
        # Use GUROBI to optimize
        model.optimize(mipsol)
        # Store results
        res = None
        if model.status == GRB.Status.OPTIMAL or model.status == GRB.Status.TIME_LIMIT:
            res = {'obj_bd': model.objBound,
                   'obj_val': model.objVal,
                   'ttb': model._ttb,
                   'sol': {j: 1 if x[j].x > 0.5 else 0 for j in x}
                  }
        return res

    def set_ds_model(self, timelimit: int=0) -> tuple:
        """Build GUROBI model for dual substitution method

        Args:
            timelimit (int, optional): time limit. No limits if leave it as default value 0.

        Returns:
            tuple: model and variable x.
        """
        N, M = self.c_low.keys(), range(self.m)
        # Build model
        model = gp.Model('ds')
        # Define variables
        x = {j: model.addVar(vtype=GRB.BINARY, name='x[{}]'.format(j)) for j in N}
        u = model.addVars(self.m, vtype=GRB.CONTINUOUS, lb=0, name='u')
        v = {j: model.addVar(vtype=GRB.CONTINUOUS, lb=0, name='v[{}]'.format(j)) for j in N}
        model.update()
        # Objective function
        if self.sense == GRB.MAXIMIZE:
            model.setObjective(gp.quicksum(self.b[i] * u[i] for i in M) + gp.quicksum(v[j] for j in N) - gp.quicksum(self.c_low[j] * x[j] for j in N), GRB.MINIMIZE)
        else:
            model.setObjective(gp.quicksum(self.b[i] * u[i] for i in M) + gp.quicksum(v[j] for j in N) + gp.quicksum(self.c_upp[j] * x[j] for j in N), GRB.MINIMIZE)
        # Add constraints
        model.addConstrs((gp.quicksum(self.a[i,j] * x[j] for j in N) <= self.b[i] for i in M), name='CAP_CON')
        if self.sense == GRB.MAXIMIZE:
            model.addConstrs((gp.quicksum(self.a[i,j] * u[i] for i in M) + v[j] + (self.c_upp[j] - self.c_low[j]) * x[j] >= self.c_upp[j] for j in N), name='DUAL_CON')
        else:
            model.addConstrs((-gp.quicksum(self.a[i,j] * u[i] for i in M) - v[j] + (self.c_low[j] - self.c_upp[j]) * x[j] <= self.c_low[j] for j in N), name='DUAL_CON')
        # Parameter setting
        model.Params.outputFlag = False
        model.Params.threads = 1
        model.Params.MIPGap = 0.0
        model.Params.lazyConstraints = 1
        if timelimit:
            model.Params.timeLimit = timelimit
        return (model, x)
    
    def write(self, filename: str):
        """write results to file

        Args:
            filename (str): file name to write
        """
        with open(filename, 'w+') as file:
            file.write("objective value: {}\n".format(self.objval))
            file.write("time to best: {:.2f}\n".format(self.ttb))
            if self.algo == 'bc':
                file.write("lower bound: {}\n".format(self.lb))
            if self.algo[:3] == 'ids':
                file.write("#iteration: {}\n".format(self.num_ite))
                file.write("best iteration: {}\n".format(self.bes_ite))
            file.write("\nSolution (non-zero variable):\n")
            for j,v in self.sol.items():
                if v < 0.5:
                    file.write(" {}\n".format(j))