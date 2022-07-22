#!/usr/bin/python
import numpy as np
import time
import math
import random
import instance as a
import logging
import networkx as nx
import util as u
import cplex
from cplex.exceptions import CplexError
import pftimeseries as pf

try:
    import gurobipy as gbp
except ImportError:
    logging.warning("Gurobi not available!!")


class ckp_sol(object):
    def __init__(self): self.obj = 0; self.solution = []; self.idx = []; self.running_time = 0; cr = None  # aproximation ratio


def OPT_v_OFFLINE(ins):
    t1 = time.time()
    m = gbp.Model("qcp")
    timehorizon = ins.timehorizon
    u.gurobi_setting(m)
    P = [0]*timehorizon
    Q = [0]*timehorizon
    x = [0] * ins.n
    v = ins.V_
    vmin = ins.v_min
    for k in ins.I: x[k] = m.addVar(vtype=gbp.GRB.BINARY, name="x_%d" % k)
    for k in ins.F: x[k] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="x_%d" % k, ub=1, lb=0)
    for turn in range(timehorizon):
        P[turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="P_%d" % turn, lb=0)
        Q[turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="Q_%d" % turn, lb=0)
    obj = gbp.quicksum(x[k] * ins.loads_utilities[k] for k in range(ins.n))
    m.setObjective(obj, gbp.GRB.MAXIMIZE)

    for turn in range(timehorizon):
        tep = gbp.LinExpr()
        tep2 = gbp.LinExpr()
        tep3 = gbp.LinExpr()
        for i in range(ins.n):
            if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]:
                cust_path = ins.customer_path_nodes[i][:-1]
                constrsum = sum([ins.topology[node][node+1]['z'][0]*ins.loads_P[i] + ins.topology[node][node+1]['z'][1]*ins.loads_Q[i] for node in cust_path])
                tep.addTerms(constrsum, x[i])
                if i in ins.F:
                    tep2.addTerms(ins.loads_P[i], x[i])
                    tep3.addTerms(1.0, x[i])
        m.addConstr(tep <= v, name="voltage_%d" % turn) #voltage constraint
        m.addConstr(tep2 <= ins.topology.graph["EVCS_capacity"], name="EV_charging_%d" % turn) #charging station capacity constraint
        m.addConstr(tep3 <= ins.topology.graph["EVSE_number"]*0.3, name="EV_charging_n_%d" % turn) #charging station slot constraint
        """
        p_sum = gbp.quicksum([x[i] * ins.loads_P[i] for i in range(ins.n) if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
        q_sum = gbp.quicksum([x[i] * ins.loads_Q[i] for i in range(ins.n) if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
        m.addQConstr(P[turn], gbp.GRB.EQUAL, p_sum)  # power flow equations
        m.addQConstr(Q[turn], gbp.GRB.EQUAL, q_sum)  # power flow equations
        m.addQConstr(P[turn]*P[turn] + Q[turn]*Q[turn], gbp.GRB.LESS_EQUAL, 74.637657*vmin, name="Current_%d" % turn)
        """
        #
    m.update()
    m.write("out.lp")

    sol = ckp_sol()
    solutionStatus = "No Solution"
    try:
        m.optimize()
        u.gurobi_handle_errors(m)
        solutionStatus = m
    except:
        print "oho"
        pass
    """
    try:
        p = cplex.Cplex()
        #p.set_log_stream(None)
        #p.set_error_stream(None)
        #p.set_warning_stream(None)
        #p.set_results_stream(None)
        p.read("out.lp")
        p.parameters.mip.display.set(1)
        #p.parameters.mip.limits.nodes.set(200000)
        p.parameters.timelimit.set(200)
        #p.parameters.mip.tolerances.integrality.set(0.01)
        #p.parameters.mip.tolerances.mipgap.set(0.2)
        #p.parameters.mip.tolerances.absmipgap.set(0.01)
        # p.parameters.mip.limits.solutions.set(1)
        p.solve()
        solutionStatus = p.solution.status[p.solution.get_status()]
        print solutionStatus
        obj = float(p.solution.get_objective_value())
    except CplexError, exc:
        print "oho"
        pass
    """

    if solutionStatus != "No Solution":
        sol.running_time = time.time() - t1
        sol.obj = obj.getValue()
        sol.idx = [k for k in range(ins.n) if x[k].x > 0]
        sol.solution = [x[k].x for k in range(ins.n)]

    return sol


def OPT_OFFLINE(ins):
    t1 = time.time()

    T = ins.topology
    m = gbp.Model("qcp")
    timehorizon = ins.timehorizon
    u.gurobi_setting(m)

    x = [0] * ins.n
    w = [0] * len(ins.F)
    v = {i: [0]*timehorizon for i in T.nodes()}
    v[0] = [ins.v_0]*timehorizon
    l = {e: [0]*timehorizon for e in T.edges()}
    P = {e: [0]*timehorizon for e in T.edges()}
    Q = {e: [0]*timehorizon for e in T.edges()}
    for k in ins.I: x[k] = m.addVar(vtype=gbp.GRB.BINARY, name="x_%d" % k)
    for k in ins.F: x[k] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="x_%d" % k, ub=1, lb=0)
    for k in range(len(ins.F)): w[k] = m.addVar(vtype=gbp.GRB.BINARY, name="w_%d" % k)
    for i in list(T.nodes())[1:]:
        for turn in range(timehorizon):
            v[i][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="v_%d_%d" % (i, turn), lb=ins.v_min, ub=ins.v_max)
    for e in T.edges():
        for turn in range(timehorizon):
            l[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="l_%s_%d" % (str(e), turn), lb=0, ub=T[e[0]][e[1]]['l']**2)
            P[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="P_%s_%d" % (str(e), turn), lb=0)
            Q[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="Q_%s_%d" % (str(e), turn), lb=0)

    obj = gbp.quicksum(x[k] * ins.loads_utilities[k] for k in range(ins.n))
    m.setObjective(obj, gbp.GRB.MAXIMIZE)
    m.update()
    for i in range(len(ins.F)):
        m.addConstr(w[i], gbp.GRB.GREATER_EQUAL, x[ins.F[i]], name="lot_%d" % i)
        m.addConstr(w[i] - 1, gbp.GRB.LESS_EQUAL, x[ins.F[i]] - 0.01, name="lotf_%d" % i)
    for turn in range(timehorizon):
        tep2 = gbp.LinExpr()
        tep = gbp.LinExpr()
        for i in range(len(ins.F)):
            if ins.loads_arrivalTime[ins.F[i]] <= turn < ins.loads_arrivalTime[ins.F[i]] + ins.loads_durations[ins.F[i]]:
                    tep2.addTerms(ins.loads_P[ins.F[i]], x[ins.F[i]])
                    tep.addTerms(1, w[i])
        m.addConstr(tep2, gbp.GRB.LESS_EQUAL, ins.topology.graph["EVCS_capacity"], name="EV_charging_%d" % turn) #charging station capacity constraint
        m.addConstr(tep, gbp.GRB.LESS_EQUAL, ins.topology.graph["EVSE_number"], name="EV_charging_n_%d" % turn) #charging station slot constraint
        for e in T.edges():
            z = T[e[0]][e[1]]['z']
            m.addQConstr(l[e][turn] * v[e[0]][turn], gbp.GRB.GREATER_EQUAL, P[e][turn] * P[e][turn] + Q[e][turn] * Q[e][turn], name="ldef_%s_%d" % (str(e), turn))  # l= |S|^2/ v_i

            psum = gbp.quicksum([x[i] * ins.loads_P[i] for i in T.node[e[1]]['N'] if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
            rhs_P = l[e][turn] * z[0] + psum + gbp.quicksum([P[h][turn] for h in list(T.edges(e[1])) if e[1] < h[1]])
            qsum = gbp.quicksum([x[i] * ins.loads_Q[i] for i in T.node[e[1]]['N'] if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
            rhs_Q = l[e][turn] * z[1] + qsum + gbp.quicksum([Q[h][turn] for h in list(T.edges(e[1])) if e[1] < h[1]])
            m.addQConstr(P[e][turn], gbp.GRB.EQUAL, rhs_P, name="P_%s_%d" % (str(e), turn))  # power flow equations
            m.addQConstr(Q[e][turn], gbp.GRB.EQUAL, rhs_Q, name="Q_%s_%d" % (str(e), turn))  # power flow equations

            rhs_v = v[e[0]][turn] + (z[0] ** 2 + z[1] ** 2) * l[e][turn] - 2 * (z[0] * P[e][turn] + z[1] * Q[e][turn])
            m.addConstr(v[e[1]][turn], gbp.GRB.EQUAL, rhs_v, name="v_%s_%d" % (str(e), turn)) # power flow equations

            if e[0] == 0:
                trafo = P[e][turn] * P[e][turn] + Q[e][turn] * Q[e][turn]
                m.addQConstr(trafo, gbp.GRB.LESS_EQUAL, T.graph["trafo"] ** 2, name="trafo_%d" % turn)  # transformer capacity

    m.update()
    m.write("out1.lp")
    sol = ckp_sol()
    solutionStatus = "No Solution"
    try:
        m.optimize()
        u.gurobi_handle_errors(m)
        solutionStatus = m
    except:
        print "oho"
        pass
    if solutionStatus != "No Solution":
        sol.running_time = time.time() - t1
        sol.obj = obj.getValue()
        sol.idx = [k for k in range(ins.n) if x[k].x > 0]
        sol.solution = [x[k].x for k in range(ins.n)]

    return sol

    """
    m.optimize()
    u.gurobi_handle_errors(m)
    sol = a.maxOPF_sol()
    sol.obj = obj.getValue()
    sol.x = {k: x[k].x for k in range(ins.n)}
    sol.idx = [k for k in ins.I if x[k].x == 1]
    sol.status = m
    

    sol = ckp_sol()
    solutionStatus = "No Solution"
    try:
        p = cplex.Cplex()
        p.set_log_stream(None)
        p.set_error_stream(None)
        p.set_warning_stream(None)
        p.set_results_stream(None)
        p.read("out.lp")
        p.parameters.mip.display.set(0)
        p.parameters.mip.limits.nodes.set(200000)
        p.parameters.timelimit.set(200)
        p.parameters.mip.tolerances.integrality.set(0.01)
        # p.parameters.mip.tolerances.mipgap.set(0.2)
        p.parameters.mip.tolerances.absmipgap.set(0.01)
        # p.parameters.mip.limits.solutions.set(1)
        p.solve()
        solutionStatus = p.solution.status[p.solution.get_status()]
        obj = float(p.solution.get_objective_value())
    except CplexError, exc:
        # print exc
        pass

    if solutionStatus != "No Solution":
        sol.running_time = time.time() - t1
        sol.obj = obj
        sol.idx = []
        for idx, val in enumerate(p.solution.get_values()):
            if val == 1:
                sol.idx.append(idx)

    return sol
    """

def OPF_ONLINE(ins, x, obj='loss', alg='FCFS'):
    T = ins.topology
    m = gbp.Model("qcp")
    timehorizon = ins.timehorizon
    u.gurobi_setting(m)

    v = {i: [0]*timehorizon for i in T.nodes()}
    v[0] = [ins.v_0]*timehorizon
    l = {e: [0]*timehorizon for e in T.edges()}
    P = {e: [0]*timehorizon for e in T.edges()}
    Q = {e: [0]*timehorizon for e in T.edges()}

    for i in list(T.nodes())[1:]:
        for turn in range(timehorizon):
            v[i][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="v_%d_%d" % (i, turn), lb=ins.v_min, ub=ins.v_max)
    for e in T.edges():
        for turn in range(timehorizon):
            l[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="l_%s_%d" % (str(e), turn), lb=0)
            P[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="P_%s_%d" % (str(e), turn), lb=0)
            Q[e][turn] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="Q_%s_%d" % (str(e), turn), lb=0)

    if obj == 'loss':
        obj = gbp.quicksum((T[0][1]['z'][0] ** 2 + T[0][1]['z'][1] ** 2) * l[(0,1)][turn] for turn in range(timehorizon))
    else:
        obj = gbp.quicksum([1])
    m.setObjective(obj, gbp.GRB.MINIMIZE)

    for turn in range(timehorizon):
        for e in T.edges():
            z = T[e[0]][e[1]]['z']
            m.addQConstr(l[e][turn] * v[e[0]][turn], gbp.GRB.GREATER_EQUAL, (P[e][turn] * P[e][turn] + Q[e][turn] * Q[e][turn]))  # l= |S|^2/ v_i

            rhs_P = l[e][turn] * z[0] + gbp.quicksum([x[i] * ins.loads_P[i] for i in T.node[e[1]]['N'] if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]]) + gbp.quicksum(
                [P[(h)][turn] for h in list(T.edges(e[1])) if e[1] < h[1]])
            rhs_Q = l[e][turn] * z[1] + gbp.quicksum([x[i] * ins.loads_Q[i] for i in T.node[e[1]]['N'] if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]]) + gbp.quicksum(
                [Q[h][turn] for h in list(T.edges(e[1])) if e[1] < h[1]])
            m.addQConstr(P[e][turn], gbp.GRB.EQUAL, rhs_P)  # power flow equations
            m.addQConstr(Q[e][turn], gbp.GRB.EQUAL, rhs_Q)  # power flow equations

            rhs_v = v[e[0]][turn] + (z[0] ** 2 + z[1] ** 2) * l[e][turn] - 2 * (z[0] * P[e][turn] + z[1] * Q[e][turn])
            m.addConstr(v[e[1]][turn], gbp.GRB.EQUAL, rhs_v) # power flow equations

            if alg == 'FCFS':
                if e[0] == 0:
                    m.addQConstr(P[e][turn] * P[e][turn] + Q[e][turn] * Q[e][turn], gbp.GRB.LESS_EQUAL, ins.capacity_profile[turn] ** 2, name="Cap_%d" % turn)  # capacity constraint

    m.update()
    m.write("out.lp")
    solutionStatus = "No Solution"
    try:
        p = cplex.Cplex()
        p.set_log_stream(None)
        p.set_error_stream(None)
        p.set_warning_stream(None)
        p.set_results_stream(None)
        p.read("out.lp")
        p.parameters.mip.display.set(0)
        p.parameters.mip.limits.nodes.set(200000)
        p.parameters.timelimit.set(200)
        p.parameters.mip.tolerances.integrality.set(0.01)
        # p.parameters.mip.tolerances.mipgap.set(0.2)
        p.parameters.mip.tolerances.absmipgap.set(0.01)
        # p.parameters.mip.limits.solutions.set(1)
        p.solve()
        solutionStatus = p.solution.status[p.solution.get_status()]
        obj = float(p.solution.get_objective_value())
    except CplexError, exc:
        # print exc
        pass

    if (solutionStatus != "No Solution") and (solutionStatus != "infeasible"):
        return True
    else:
        return False


def min_loss_OPF(ins, x, cons=''):
    t1 = time.time()
    T = ins.topology
    m = gbp.Model("qcp")
    u.gurobi_setting(m)

    dummy_p = {i: 0 for i in T.nodes()}
    dummy_q = {i: 0 for i in T.nodes()}
    v = {i: 0 for i in T.nodes()}
    v[0] = ins.v_0
    l = {e: 0 for e in T.edges()}
    P = {e: 0 for e in T.edges()}
    Q = {e: 0 for e in T.edges()}
    for k in T.nodes(): dummy_p[k] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="dummy_p[%d]" % k, lb=0.005, ub=0.006)
    for k in T.nodes(): dummy_q[k] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="dummy_q[%d]" % k, lb=0.005, ub=0.006)
    for k in T.nodes()[1:]:
        v[k] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="v_%d" % k)
    for e in T.edges(): l[e] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="l_%s" % str(e))
    for e in T.edges(): P[e] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="P_%s" % str(e))
    for e in T.edges(): Q[e] = m.addVar(vtype=gbp.GRB.CONTINUOUS, name="Q_%s" % str(e))
    m.update()

    #obj = gbp.quicksum([l[e] * np.sqrt(T[e[0]][e[1]]['z'][0] ** 2 + T[e[0]][e[1]]['z'][0] ** 2) for e in T.edges()])
    obj = gbp.quicksum([1])
    m.setObjective(obj, gbp.GRB.MINIMIZE)
    pesum = gbp.LinExpr()
    qesum = gbp.LinExpr()
    for e in T.edges():
        z = T[e[0]][e[1]]['z']

        m.addQConstr(l[e] * v[e[0]], gbp.GRB.GREATER_EQUAL, (P[e] * P[e] + Q[e] * Q[e]),
                     "l_%s" % str(e))  # l= |S|^2/ v_i
        #m.addConstr(dummy_p[e[1]] >= 0, "dummy_P_%d" % e[1])
        # index_set = set(T.node[e[1]]['N']).intersection(set(idx))
        rhs_P = l[e] * z[0] + sum([ins.loads_P[k] * x[k] for k in T.node[e[1]]['N']]) + gbp.quicksum(
            [P[(e[1], h)] for h in T.edge[e[1]].keys() if e[1] < h]) #+ dummy_p[e[1]]
        rhs_Q = l[e] * z[1] + sum([ins.loads_Q[k] * x[k] for k in T.node[e[1]]['N']]) + gbp.quicksum(
            [Q[(e[1], h)] for h in T.edge[e[1]].keys() if e[1] < h]) #+ dummy_q[e[1]]
        m.addQConstr(P[e], gbp.GRB.EQUAL, rhs_P, "P_%s=" % str(e))
        m.addQConstr(Q[e], gbp.GRB.EQUAL, rhs_Q, "Q_%s=" % str(e))

        rhs_v = v[e[0]] + (z[0] ** 2 + z[1] ** 2) * l[e] - 2 * (z[0] * P[e] + z[1] * Q[e])
        m.addConstr(v[e[1]], gbp.GRB.EQUAL, rhs_v, "v_%d=" % e[1])
        pesum += P[e]
        qesum += Q[e]
        #if cons == 'C' or cons == '':
            #if T[e[0]][e[1]]['C'] == 4:
                #m.addQConstr(P[e] * P[e] + Q[e] * Q[e], gbp.GRB.LESS_EQUAL, T[e[0]][e[1]]['C'] ** 2,
                             #"C_%s" % str(e))  # capacity constraint
        if cons == 'V' or cons == '':
            m.addConstr(v[e[1]], gbp.GRB.GREATER_EQUAL, ins.v_min, "v_%d" % e[1])  # voltage constraint

    #m.addQConstr(pesum * pesum + qesum * qesum, gbp.GRB.LESS_EQUAL, T[0][1]['C'] ** 2, "C_%s" % str(e))  # capacity constraint
    m.update()
    m.write("out1.lp")
    m.optimize()
    sol = a.maxOPF_sol()
    sol.status = m.status
    sol.m = m
    sol.running_time = time.time() - t1

    if (u.gurobi_handle_errors(m) == False):
        sol.obj = -np.inf
        return None
    sol.obj = obj.getValue()
    for i in m.getVars():
        if i.VarName == "v_4":
            return i.X


def OPT_gurobi_online(ins, C, timehorizon):
    t1 = time.time()
    MGT = ins.topology
    v = {i: 0 for i in MGT.nodes()}
    v[0] = ins.v_0
    v_min = ins.v_min
    m = gbp.Model("qcp")
    m.setParam("OutputFlag", 0)
    m.setParam("TimeLimit", 200)
    # m.setParam("MIPGapAbs", 0.000001)
    m.setParam("MIPGapAbs", 0)
    # m.setParam("MIPGap", 0)
    # m.setParam("SolutionLimit", 1)
    m.setParam("IntFeasTol", 0.00001)
    C = np.array(C).astype(dtype='int64')
    x = [None] * ins.n
    for i in range(ins.n):
        x[i] = m.addVar(vtype=gbp.GRB.BINARY, name="x_%d" % i)
    m.update()

    obj = gbp.quicksum(x[c] * ins.loads_utilities[c] for c in range(ins.n))
    m.setObjective(obj, gbp.GRB.MAXIMIZE)

    for turn in range(timehorizon):
        if [ins.loads_S[i] for i in range(ins.n) if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]]:
            #constrlist = []
            tep = gbp.LinExpr()
            for i in range(ins.n):
                constrsum = 0
                if ins.loads_arrivalTime[i] <= turn < ins.loads_arrivalTime[i] + ins.loads_durations[i]:
                    constrsum += sum([ins.topology[node][node+1]['z'][0]*ins.loads_P[i] + ins.topology[node][node+1]['z'][1]*ins.loads_Q[i] for node in ins.customer_path_nodes[i] if node < 4])
                #constrlist.append(x[i] * constrsum)
                    tep.addTerms(constrsum, x[i])
            #constraint = gbp.quicksum(constrlist)
            m.addConstr(tep <= 0.5*(v[0] - v_min))
    m.update()
    m.write("model.lp")
    sol = ckp_sol()
    solutionStatus = "No Solution"
    try:
        p = cplex.Cplex()
        p.read("model.lp")
        p.parameters.mip.limits.nodes.set(200000)
        p.parameters.timelimit.set(200)
        p.parameters.mip.tolerances.integrality.set(0.01)
        #p.parameters.mip.tolerances.mipgap.set(0.2)
        p.parameters.mip.tolerances.absmipgap.set(0.01)
        #p.parameters.mip.limits.solutions.set(1)
        p.parameters.mip.display.set(0)
        p.solve()
        solutionStatus = p.solution.status[p.solution.get_status()]
        obj = float(p.solution.get_objective_value())
    except CplexError, exc:
        #print exc
        pass
    if solutionStatus != "No Solution":
    #m.optimize()
    #if m.status != gbp.GRB.status.OPTIMAL:
        #print "!!!!!!! Gurobi returned non-OPT solution !!!!!!!!"
        #print " status: ", m.status
        sol.running_time = time.time() - t1
        sol.obj = obj#obj.getValue()
        sol.idx = []

        for idx, val in enumerate(p.solution.get_values()):
            if val == 1:
                sol.idx.append(idx)
        #assert sum([ins.loads_utilities[idx] for idx in sol.idx]) == sol.obj
        return sol


def PD(timehorizon, capacities, theset, coefficients, intervals, x, y, index, utility, type, curround):
    amax = [0]*timehorizon
    for t in intervals[curround]:
        amax[t] = max([coefficients[j][t] for j in range(len(theset))])
    tmax = max([len(intervals[k]) for k in theset])
    while sum([y[t]*coefficients[index][t] for t in intervals[curround]]) < 1:
        if type == "s":
            x[index] += 0.4*utility
        else:
            x[index] += 0.05*utility
        for t in intervals[curround]:
            b = (1./(tmax*amax[t])) * (math.exp((1./(2*capacities[t]))*sum([x[j] * coefficients[j][t] for j in range(len(theset))]))-1)
            y[t] = max(y[t], b)
    return x[index], y


def PDAOnline(ins, timehorizon, alpha, delta, PPMG):
    initialtimehorizon = timehorizon
    t1 = time.time()
    #Cprime = np.array([i*ins.topology.graph['S_base'] for i in ins.capacity_profile]).astype(dtype='int64')
    y = [0] * (timehorizon+ins.n)
    x = [0] * ins.n
    rx = [0] * ins.n
    ytilde = [0] * (timehorizon)
    xtilde = [0] * ins.n
    xhat = [0] * ins.n
    xf = [0] * len(ins.F)
    yf = [0] * (2*timehorizon + len(ins.F))
    Is = []
    Il = []
    F = []
    s = 0
    l = 0
    f = 0
    T = []
    a = [[0 for i in range(timehorizon+ins.n)] for j in range(ins.n)]
    u = [[0 for i in range(timehorizon)] for j in range(ins.n)]
    #C = [int(i*ins.topology.graph['S_base']) for i in ins.capacity_profile]
    b = [ins.V_/2.]*timehorizon
    bf = np.append(b, [ins.topology.graph["EVCS_capacity"]]*timehorizon)
    af = [[0 for i in range(2*timehorizon+len(ins.F))] for j in range(len(ins.F))]
    ftimehorizon = 2*timehorizon

    for curround in range(ins.n):
        T.append(list(range(int(ins.loads_arrivalTime[curround]), int(ins.loads_arrivalTime[curround] + ins.loads_durations[curround]))))
        if curround in ins.F:
            #fractional demands
            F.append(curround)
            T[curround].extend(list(range(int(initialtimehorizon + ins.loads_arrivalTime[curround]), int(initialtimehorizon + ins.loads_arrivalTime[curround] + ins.loads_durations[curround]))))
            T[curround].append(ftimehorizon)
            bf = np.append(bf, [ins.loads_utilities[curround]])
            for t in T[curround]:
                if t < initialtimehorizon:
                    af[f][t] = ins.loads_G[curround] / (1. * ins.loads_utilities[curround])
                else:
                    af[f][t] = ins.loads_S[curround] / (1. * ins.loads_utilities[curround])
            af[f][ftimehorizon] = 1
            ftimehorizon += 1
            xf[f], yf = PD(ftimehorizon, [bf[t] for t in range(ftimehorizon)], F, af, T, xf, yf, f, ins.loads_utilities[curround], "f", curround)
            f += 1
            temp = np.array(af)
            themin = np.amin(temp[temp[:] != 0])
            f_scalefactor = 2. * math.log(1 + ftimehorizon * np.amax(af) / (1.0 * themin), 2) * ins.loads_utilities[curround]
            rx[curround] = xf[f - 1] / f_scalefactor
        else:
            if ins.loads_G[curround] <= delta*min(b[:initialtimehorizon]):
                #small demands
                Is.append(curround)
                T[curround].append(timehorizon)
                b = np.append(b, [2*ins.loads_utilities[curround]])
                for t in T[curround]:
                    a[s][t] = ins.loads_G[curround]/(1.*ins.loads_utilities[curround])
                a[s][timehorizon] = 1
                timehorizon += 1
                xhat[s], y = PD(timehorizon, [b[t]/2. for t in range(timehorizon)], Is, a, T, xhat, y, s, ins.loads_utilities[curround], "s", curround)
                #xhat[s], y = PD(timehorizon, [C[t] for t in range(timehorizon)], Is, a, T, xhat, y, s, ins.loads_utilities[curround], "s")
                s += 1
            else:
                #largr demands
                Il.append(curround)
                for t in T[curround]:
                    u[l][t] = 1./ins.loads_utilities[curround]
                xtilde[l], ytilde = PD(initialtimehorizon, [1]*initialtimehorizon, Il, u, T, xtilde, ytilde, l, ins.loads_utilities[curround], "l", curround)
                l += 1

            temp = np.array(a)
            themin = np.amin(temp[temp[:] != 0])
            s_scalefactor = 2. * math.log(1 + timehorizon * np.amax(a) / (1.0 * themin), 2) * ins.loads_utilities[curround]
            l_scalefactor = 2. * math.log(1 + initialtimehorizon * max(ins.loads_utilities[:curround+1]) / (1.0 * min(ins.loads_utilities[:curround+1])), 2) * ins.loads_utilities[curround]
            if curround in Is:
                x[curround] = xhat[s-1]/s_scalefactor
                roundprb = random.uniform(0.09, 0.8)
            else:
                x[curround] = xtilde[l - 1]/l_scalefactor
                roundprb = random.uniform(0, 0.00001)

            if roundprb <= x[curround]:
                rx[curround] = 1

        #feasibility check
        if rx[curround] > 0:
            minvoltage, max_line_loading, max_trafo_loading = pf.simulate_MG(PPnet=PPMG, topology=ins.topology, timesteps=initialtimehorizon,
                                                                             sol=rx, ins=ins, mode="subroutine",
                                                                             curround=curround,
                                                                             turn=ins.loads_arrivalTime[curround])
            if minvoltage < 0.95 or max_line_loading > 100 or max_trafo_loading > 100:
                rx[curround] = 0
            else:
                term2 = sum([rx[i] * ins.loads_S[i] for i in range(curround + 1) if
                             i in ins.F and ins.loads_arrivalTime[i] <= ins.loads_arrivalTime[curround] <
                             ins.loads_arrivalTime[i] + ins.loads_durations[i]])
                if term2 > ins.topology.graph["EVCS_capacity"]:
                    rx[curround] = 0
                    continue
                term3 = sum([1 for i in range(curround + 1) if i in ins.F and rx[i] > 0 and ins.loads_arrivalTime[i] <= ins.loads_arrivalTime[curround] < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
                if term3 > ins.topology.graph["EVSE_number"]:
                    rx[curround] = 0
                    continue

    sol = ckp_sol()
    sol.running_time = time.time() - t1
    sol.obj = sum([rx[i]*ins.loads_utilities[i] for i in range(ins.n)])
    sol.idx = [k for k in range(ins.n) if rx[k] > 0]
    sol.solution = rx
    return sol


def FirstComefirstserve(ins=None, PPMG=None, MGTree=None, T=None):

    t1 = time.time()
    x = [0] * ins.n
    for curround in range(ins.n):
        x[curround] = 1
        minvoltage, max_line_loading, max_trafo_loading = pf.simulate_MG(PPnet=PPMG, topology=MGTree, timesteps=T, sol=x, ins=ins, mode="subroutine", curround=curround, turn=ins.loads_arrivalTime[curround])
        if minvoltage < 0.95 or max_line_loading > 100 or max_trafo_loading > 100:
            x[curround] = 0
        else:
            term2 = sum([x[i] * ins.loads_S[i] for i in range(curround + 1) if i in ins.F and ins.loads_arrivalTime[i] <= ins.loads_arrivalTime[curround] < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
            if term2 > MGTree.graph["EVCS_capacity"]:
                x[curround] = 0
                continue
            term3 = sum([x[i] for i in range(curround + 1) if i in ins.F and ins.loads_arrivalTime[i] <= ins.loads_arrivalTime[curround] < ins.loads_arrivalTime[i] + ins.loads_durations[i]])
            if term3 > MGTree.graph["EVSE_number"]:
                x[curround] = 0

    """
    t1 = time.time()
    x = [0]*ins.n
    v = ins.V_
    for curround in range(ins.n):
        cust_path = ins.customer_path_nodes[curround][:-1]
        tmp = sum([ins.topology[node][node + 1]['z'][0] * ins.loads_P[curround] + ins.topology[node][node + 1]['z'][1] *
                    ins.loads_Q[curround] for node in cust_path])
        for turn in range(curround):
            if ins.loads_arrivalTime[turn] <= curround < ins.loads_arrivalTime[turn] + ins.loads_durations[turn]:
                if x[turn] != 0:
                    cust_path = ins.customer_path_nodes[turn][:-1]
                    tmp += sum([ins.topology[node][node+1]['z'][0]*ins.loads_P[turn] + ins.topology[node][node+1]['z'][1]*ins.loads_Q[turn] for node in cust_path])
                    print tmp
        print tmp, v
        if curround == 5:
            exit()
        if tmp <= v:
            x[curround] = 1
            """

    sol = ckp_sol()
    sol.running_time = time.time() - t1
    sol.obj = sum([x[i]*ins.loads_utilities[i] for i in range(ins.n)])
    sol.idx = [k for k in range(ins.n) if x[k] == 1]
    sol.solution = x
    return sol


