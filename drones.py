from gekko import GEKKO
import random
import timeit


m = GEKKO(remote=True)
times = []
for n in range(10, 100, 10):
    start = timeit.default_timer()
    x = [[m.Var(0, lb = 0, ub = 1, integer = True) for i in range(n)] for j in range(n)]
    g = [[m.Var(0, lb = 0, ub = n, integer = True) for i in range(n)] for j in range(n)]
    c = [[random.randint(1, 10 * n) for i in range(n)] for j in range(n)]
    for i in range(n):
        m.Equation(sum([x[i][j] for j in range(n) if i != j]) == 1)
        m.Equation(sum([x[j][i] for j in range(n) if i != j]) == 1)
    for i in range(1, n):
        m.Equation(sum([g[j][i] for j in range(n)]) - sum([g[i][j] for j in range(1, n)]) == 1)

    for i in range(n):
        for j in range(1, n):
            m.Equation(g[i][j] <= (n - 1) * x[i][j])

    m.Obj(sum([sum([c[i][j] * x[i][j] for j in range(n) if j != i]) for i in range(n)]))
    m.options.SOLVER=1
    m.solve(disp=True)
    # print(x)
    # print(c)
    stop = timeit.default_timer()
    print('Time for ' + str(n) + ' vertices: ', stop - start)
    times.append(stop - start)

print(times)

# print(x)
# print(u)

# u = 0, 1, 2
# 1 -> 0,1 1,2 2,3