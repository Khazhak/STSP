from gurobipy import Model, GRB, quicksum
import numpy as np


class TSPSolver():
    def __init__(self, n, c, v_lower = None):
        self.n = n
        self.c = c
        self.v_lower = v_lower

    def solve(self, algorithm = 'P1', N_lower = None, v_vals = None):
        if algorithm == 'P2':
            m = Model()
            m.Params.LogToConsole = 0
            N = [i for i in range(self.n)]
            N_ = [i for i in range(1, self.n)]
            A = [(i, j) for i in N for j in N]
            x = [[0 for i in N] for j in N]
            g = [[0 for i in N] for j in N]


            for (i, j) in A:
                x[i][j] = m.addVar(vtype=GRB.BINARY, name="x_%d_%d" % (i, j))
            for (i, j) in A:
                g[i][j] = m.addVar(vtype=GRB.CONTINUOUS, name="g_%d_%d" % (i, j))


            for i in N:
                m.addConstr(x[i][i] == 0)
                m.addConstr(quicksum(x[i][j] for j in N if j != i) == 1)
                m.addConstr(quicksum(x[j][i] for j in N if j != i) == 1)

            for i in N:
                if i != 0:
                    m.addConstr(quicksum(g[j][i] for j in N) - quicksum(g[i][j] for j in N_) == 1)
                for j in N_:
                    if i != j:
                        m.addConstr(g[i][j], GRB.LESS_EQUAL, (self.n - 1) * x[i][j])
                        m.addConstr(0 <= g[i][j])
            

            m.setObjective(quicksum(x[i][j] * self.c[i][j] for i, j in A), 0)


            m.modelSense = GRB.MINIMIZE
            m.optimize()
            x_array = [[x[i][j].x for i in range(self.n)] for j in range(self.n)] 

            return x_array, m.objVal
        elif algorithm == 'P1':
            if N_lower is None:
                print('N_lower should be given if algorithm=P1')
                raise NotImplementedError

            m = Model()
            m.Params.LogToConsole = 0
            m.Params.NonConvex = 2
            N = [i for i in range(self.n)]
            N_ = [i for i in range(1, self.n)]
            A = [(i, j) for i in N for j in N]
            x = [[0 for i in N] for j in N]
            g = [[0 for i in N] for j in N]
            v = [[0 for i in N] for j in N]
            v_div = [[0 for i in N] for j in N]

            for (i, j) in A:
                x[i][j] = m.addVar(vtype=GRB.BINARY, name="x_%d_%d" % (i, j))
            for (i, j) in A:
                g[i][j] = m.addVar(vtype=GRB.CONTINUOUS, name="g_%d_%d" % (i, j))
            for (i, j) in A:
                v[i][j] = m.addVar(vtype=GRB.CONTINUOUS, name="v_%d_%d" % (i, j), lb = self.v_lower[i][j], ub = 1)
            for (i, j) in A:
                v_div[i][j] = m.addVar(vtype=GRB.CONTINUOUS, name="v_div_%d_%d" % (i, j), ub = 1 / self.v_lower[i][j], lb = 1.0)


            for i, j in A:
                m.addConstr(v[i][j] * v_div[i][j] == 1)

            for i in N:
                m.addConstr(x[i][i] == 0)
                m.addConstr(quicksum(x[i][j] for j in N if j != i) == 1)
                m.addConstr(quicksum(x[j][i] for j in N if j != i) == 1)

            for i in N:
                if i != 0:
                    m.addConstr(quicksum(g[j][i] for j in N) - quicksum(g[i][j] for j in N_) == 1)
                for j in N_:
                    if i != j:
                        m.addConstr(g[i][j], GRB.LESS_EQUAL, (self.n - 1) * x[i][j])
                        m.addConstr(0 <= g[i][j])

            m.addConstr(quicksum(x[i][j] * (1 - v[i][j]) for i, j in A) >= N_lower)

            m.setObjective(quicksum(x[i][j] * self.c[i][j] * v_div[i][j] for i, j in A), 0)

            m.modelSense = GRB.MINIMIZE
            m.optimize()
            m.write('./m.lp')

            if m.status == GRB.OPTIMAL:
                x_array = [[x[i][j].x for i in range(self.n)] for j in range(self.n)] 
                v_array = [[v[i][j].x for i in range(self.n)] for j in range(self.n)] 

                return x_array, v_array, m.objVal
            return None, None, None
        elif algorithm == 'RP1':
            if N_lower is None:
                print('N_lower should be given if algorithm=opt')
                raise NotImplementedError
            if v_vals is None:
                print('v_vals should be given if algorithm=opt')
                raise NotImplementedError

            m = Model()
            m.Params.LogToConsole = 0

            N = [i for i in range(self.n)]
            N_ = [i for i in range(1, self.n)]
            V = [i for i in  range(len(v_vals[0][0]))]
            A = [(i, j) for i in N for j in N]
            A_V = [(i, j, k) for i in N for j in N for k in V]
            x = [[0 for i in N] for j in N]
            g = [[0 for i in N] for j in N]
            v_id = [[[0 for i in V] for j in N] for k in N]


            for (i, j) in A:
                x[i][j] = m.addVar(vtype=GRB.BINARY, name="x_%d_%d" % (i, j))
            for (i, j) in A:
                g[i][j] = m.addVar(vtype=GRB.CONTINUOUS, name="g_%d_%d" % (i, j))
            for (i, j, k) in A_V:
                v_id[i][j][k] = m.addVar(vtype=GRB.BINARY, name="v_id_%d_%d_%d" % (i, j, k))



            for i in N:
                m.addConstr(x[i][i] == 0)
                m.addConstr(quicksum(x[i][j] for j in N if j != i) == 1)
                m.addConstr(quicksum(x[j][i] for j in N if j != i) == 1)
            for i in N:
                for j in N:
                    m.addConstr(quicksum(v_id[i][j][k] for k in V) == 1)


            for i in N:
                if i != 0:
                    m.addConstr(quicksum(g[j][i] for j in N) - quicksum(g[i][j] for j in N_) == 1)
                for j in N_:
                    if i != j:
                        m.addConstr(g[i][j], GRB.LESS_EQUAL, (self.n - 1) * x[i][j])
                        m.addConstr(0 <= g[i][j])

            m.addConstr(quicksum(x[i][j] * v_id[i][j][k] * (1 - v_vals[i][j][k]) for i, j, k in A_V) >= N_lower)

            m.setObjective(quicksum(x[i][j] * self.c[i][j] * v_id[i][j][k] / v_vals[i][j][k] for i, j, k in A_V), 0)

            m.modelSense = GRB.MINIMIZE
            m.optimize()

            opt_relax_N = 0
            for i in range(self.n):
                for j in range(self.n):
                    for k in range(len(v_vals[i][j])):
                        opt_relax_N += (1 - v_vals[i][j][k]) * x[i][j].x * v_id[i][j][k].x
            m.write('./m.lp')

            if m.status == GRB.OPTIMAL:
                x_array = [[x[i][j].x for i in N] for j in N] 
                v_id_array = [[[v_id[i][j][k].x for k in V] for j in N] for i in N] 

                return x_array, v_id_array, m.objVal, opt_relax_N
            return None, None, None
        else:
            raise NotImplementedError
            
