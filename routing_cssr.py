import sys
import time
import pyomo.environ as pyo
from pyomo.opt import SolverFactory
from models_cssr import load_case
from models_cssr import deterministic, subproblem, master_problem, solve_model
import numpy as np
from plot import plotingLines, plotingBuses
import os


def robust_optimization(data, gama):
    # log file
    flog = open('routing.log', 'w')
    flog.write(f'case_file: {CASE_FILE}\n')

    # define models
    DP = deterministic(data)
    SP = subproblem(data, gama)
    MP = master_problem(data)

    # SOLVE DETERMINISTIC MODEL
    opt = SolverFactory('cplex')
    t0 = time.time()
    opt.options['mip tolerances absmipgap'] = 1e-04
    opt.options['mip tolerances mipgap'] = 1e-02
    opt.options['mip tolerances integrality'] = 1e-03
    opt.options['mip display'] = 2
    opt.solve(DP, tee=False)
    t = time.time()-t0
    print(f'deterministic solved, t = {t:.2f} sec')
    flog.write(f'deterministic solved, t = {t:.2f} sec\n')
    flog.write(f'\ndeterministic, DP.obj = {pyo.value(DP.obj):.2f}\n')
    flog.write(f'deterministic, DP.Cc  = {pyo.value(DP.Cc):.2f}\n')
    flog.write(f'deterministic, DP.Cl  = {pyo.value(DP.Cl):.2f}\n')
    v_det = np.zeros((data['A'].shape[1], data['nc']))
    for i in DP.lines:
        for j in DP.cables:
            v_det[i, j] = round(pyo.value(DP.v[i, j]))
    flog.write(f'deterministic, v_det = {v_det}\n')
    P_det = np.zeros(data['nl'])
    for i in DP.lines:
        P_det[i] = pyo.value(DP.P[i])
    flog.write(f'deterministic, P_det = {P_det}\n')
    V_det = np.zeros(data['nb'])
    for i in DP.buses:
        V_det[i] = np.sqrt(pyo.value(DP.W[i]))
    flog.write(f'deterministic, V_det = {V_det}\n')

    # iterations
    it = 0
    UB, LB = 1e10, -1e10
    eps = 0.01
    while abs(UB-LB) > eps and it <= 10:
        it = it+1
        print(f'\niteration = {it}')

        if it == 1:
            # take cable selection from the solution of the deterministic model
            for i in SP.lines:
                for j in SP.cables:
                    SP.v[i, j] = pyo.value(DP.v[i, j])
        else:
            # get cable selection from the solution of the master problem and use it in the subproblem
            # print('# SP.v[i,j] = round(pyo.value(MP.v[i,j]))')
            for i in SP.lines:
                for j in SP.cables:
                    SP.v[i, j] = round(pyo.value(MP.v[i, j]))

        # SOLVE THE SUBPROBLEM
        t0 = time.time()
        solver = 'cplex'
        opt = SolverFactory(solver)
        opt.options['mip tolerances absmipgap'] = 1e-04
        opt.options['mip tolerances mipgap'] = 1e-02
        opt.options['mip tolerances integrality'] = 1e-03
        opt.options['optimalitytarget'] = 3
        opt.options['mip display'] = 2
        solve_model(opt, SP)
        UB = pyo.value(SP.Cc + SP.obj)
        t = time.time() - t0
        print(f'subproblem solved, t = {t:.2f} sec')
        flog.write(f'\niter = {it}, SP.obj = {pyo.value(SP.obj):.2f}\n')
        flog.write(f'iter = {it}, SP.Cc = {pyo.value(SP.Cc):.2f}\n')
        flog.write(f'iter = {it}, SP.Cl = {pyo.value(SP.Cl):.2f}\n')
        flog.write(f'iter = {it}, UB = {UB:.2f}\n')
        Pd = [round(pyo.value(SP.Pd[i])*1000, 2) for i in SP.loads]
        Pd_rel = [round(Pd[i-1]/(data['Pmax'][i]*1000), 2) for i in SP.loads]
        flog.write(f'iter = {it}, Pd = {Pd}\n')
        flog.write(f'iter = {it}, Pd/Pmax = {Pd_rel}\n')
        flog.write(f't = {t:.2f} sec\n')

        # add new variables and constraints (C&CG)
        if it > 1:
            # for it == 1 variables are added when MP is created
            MP.it.add(it)
            for i in MP.lines:
                MP.U.add((i, it))
            for i in MP.buses:
                MP.W.add((i, it))

        # get load demand from the subproblem and use it in the master problem
        for i in MP.loads:
            MP.Pd[i, it] = SP.Pd[i]
            MP.Qd[i, it] = SP.Qd[i]
        # calculate load flow in each branch
        A2 = data['A'][1:, :].todense()
        for i in MP.lines:
            MP.P[i, it] = -sum(np.linalg.inv(A2)[i, j] *
                               MP.Pd[j+1, it] for j in MP.lines)
            MP.Q[i, it] = -sum(np.linalg.inv(A2)[i, j] *
                               MP.Qd[j+1, it] for j in MP.lines)
            MP.S[i, it] = (MP.P[i, it]**2 + MP.Q[i, it]**2)**0.5

        # add limit on eta using optimal value from the subproblem
        constant = 8760*data['cl']*data['beta']*1000/data['Vs']**2
        MP.etaLimit.add(expr=MP.eta >= sum(sum(MP.v[i, j]*data['r'][j] for j in MP.cables) * data['d'][i]*(
            MP.P[i, it]**2+MP.Q[i, it]**2)*constant for i in MP.lines))

        # supply bus voltage
        MP.supplyBus.add(expr=MP.W[0, it] == data['Vs']**2)

        # line voltage equation
        for i in MP.lines:
            MP.lineVoltage.add(expr=MP.U[i, it] == sum(
                data['A'][j, i]*MP.W[j, it] for j in MP.buses))

        # square voltage losses
        for i in MP.lines:
            MP.voltageLosses.add(expr=MP.U[i, it] == 2*(MP.P[i, it]*sum(MP.v[i, j]*data['r'][j]
                                 for j in MP.cables)+MP.Q[i, it]*sum(MP.v[i, j]*data['x'][j] for j in MP.cables))*data['d'][i])
        # max power flow in line
        for i in MP.lines:
            MP.maxPowerFlow.add(expr=MP.S[i, it] <= sum(
                MP.v[i, j]*data['Smax'][j] for j in MP.cables))

        # SOLVE MASTER PROBLEM
        t0 = time.time()
        solver = 'cplex'
        opt = SolverFactory(solver)
        opt.options['mip tolerances absmipgap'] = 1e-04
        opt.options['mip tolerances mipgap'] = 1e-02
        opt.options['mip tolerances integrality'] = 1e-03
        opt.options['mip display'] = 2
        opt.options['optimalitytarget'] = 0
        solve_model(opt, MP)
        LB = pyo.value(MP.obj)
        t = time.time() - t0
        print(f'master problem solved, t = {t:.2f} sec')
        flog.write(f'\niter = {it}, MP.Cc = {pyo.value(MP.Cc):.2f}\n')
        flog.write(f'iter = {it}, MP.eta = {pyo.value(MP.eta):.2f}\n')
        flog.write(f'iter = {it}, MP.obj = {pyo.value(MP.obj):.2f}\n')
        flog.write(f'iter = {it}, LB = {LB:.2f}\n')
        flog.write(f't = {t:.2f} sec\n')
        v = np.zeros((data['A'].shape[1], data['nc']))
        for i in DP.lines:
            for j in DP.cables:
                v[i, j] = round(pyo.value(MP.v[i, j]))
        flog.write(f'iter = {it}, v = {v}\n')

    flog.write('\ndifference robust/deterministic\n')
    obj_diff = (LB/pyo.value(DP.obj) - 1)*100
    flog.write(f'(MP.obj/DP.obj - 1)={obj_diff: .2f} %\n')
    flog.close()

    # PLOT
    # Ploting deterministic solution
    lineSize = np.zeros(data['A'].shape[1])
    for i in range(data['A'].shape[1]):
        k = np.where(v[i, :] == 1)
        lineSize[i] = k[0]
    location = os.path.join("..", "Sliki")
    nameLines = 'graph-lines-case85-g1.txt' 
    plotingLines(data, location, nameLines, lineSize)
    nameBuses = 'graph-buses-case85.txt'
    plotingBuses(data, location, nameBuses)


if __name__ == '__main__':

    if len(sys.argv) == 1:
        CASE_FILE = 'case85_cssr.py'
    else:
        CASE_FILE = sys.argv[1]

    DATA = load_case(CASE_FILE)

    robust_optimization(DATA, 1)
