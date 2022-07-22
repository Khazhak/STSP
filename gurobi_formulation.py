from gurobipy import Model, GRB, quicksum
import random
import time
import numpy as np
from tsp_solver import TSPSolver

results = []
times = []

def dist(x1, y1, x2, y2):
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5

for n in range(10, 41, 3):
    tm = []
    x_cord = [random.random() * 10 * n for i in range(n)]
    y_cord = [random.random() * 10 * n for i in range(n)]
    c = [[dist(x_cord[i], y_cord[i], x_cord[j], y_cord[j]) for i in range(n)] for j in range(n)]
    v_lower = [[0.1 + random.random() * 0.9 for i in range(n)] for j in range(n)]
    v = [[v_lower[i][j] for i in range(n)] for j in range(n)]
    v_split = 10
    v_vals = [[np.arange(v_lower[i][j], 1.00001, (1 - v_lower[i][j])/ (v_split - 1)) for i in range(n)] for j in range(n)]
    N_lower = n / 4

    solver = TSPSolver(n, c, v_lower)

    start = time.time()

    x, _ = solver.solve(algorithm='tsp')

    tsp_objective = 0
    tsp_N = 0
    cnt = 1000

    for i in range(n):
        for j in range(n):
            # if x[i][j] > 0.5:
            #     print('tsp: ', i, j, v[i][j])
            tsp_objective += x[i][j] * c[i][j] / v[i][j]
            tsp_N += (1 - v[i][j]) * x[i][j]
    
    print('objective of the tsp algorithm at the beggining: ', tsp_objective)
    while cnt > 0 and tsp_N > N_lower + 1e-4:
        tsp_N = 0
        tsp_objective = 0
        i = random.randint(0, n - 1)
        j = random.randint(0, n - 1)
        if v[i][j] < 1:
            v[i][j] = max(1, v[i][j] + 0.1)

        for i in range(n):
            for j in range(n):
                tsp_objective += x[i][j] * c[i][j] / v[i][j]
                tsp_N += (1 - v[i][j]) * x[i][j]
        cnt -= 1
    
    print('objective of the tsp algorithm: ', tsp_objective)
    print('N of the tsp algorithm: ', tsp_N)
    end = time.time()
    tm.append(end-start)
    start = time.time()
    if n < 11:
        x, v, opt_objective = solver.solve(algorithm='opt', N_lower=N_lower)

        opt_N = 0
        for i in range(n):
            for j in range(n):
                # if x[i][j] > 0.5:
                #     print('opt: ', i, j, v[i][j])
                opt_N += (1 - v[i][j]) * x[i][j]

        print('objective of the opt algorithm: ', opt_objective)
        print('N of the opt algorithm: ', opt_N)
    else:
        opt_objective = -1
        opt_N = -1
    end = time.time()
    tm.append(end-start)
    start = time.time()

    x, v_id, opt_relax_objective, opt_relax_N_old = solver.solve(algorithm='opt_relax', N_lower=N_lower, v_vals = v_vals)

    end = time.time()
    tm.append(end-start)
    # opt_relax_N = 0
    # for i in range(n):
    #     for j in range(n):
    #         for k in range(len(v_vals)):
    #             if x[i][j] > 0.5 and v_id[i][j][k] > 0.5:
    #                 print('opt_relax: ', i, j, k, v_vals[k])
    #             opt_relax_N += (1 - v_vals[i][j][k]) * x[i][j] * v_id[i][j][k]
    

    print('objective of the opt_relax algorithm: ', opt_relax_objective)
    # print('N of the opt_relax algorithm: ', opt_relax_N)
    print('N of the opt_relax old algorithm: ', opt_relax_N_old)
    results.append([tsp_objective, tsp_N, opt_objective, opt_N, opt_relax_objective, opt_relax_N_old])
    times.append(tm)
    print(tm)

np.save( './results', np.array(results))
np.save( './times', np.array(times))
