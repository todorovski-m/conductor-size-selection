from email.policy import default
import importlib
from msilib.schema import Binary
from tkinter import END
import numpy as np
import networkx as nx
import pyomo.environ as pyo
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt


def solve_model(o, m):
    if hasattr(o, 'set_instance'):
        o.set_instance(m)  # required by gurobi
    r = o.solve(m, tee=False)
    return r


def load_case(case_file):
    # reading input file
    print(f'\ncase_file={case_file}')
    case_file = case_file.replace('.py', '')
    my_module = importlib.import_module(case_file)
    ds = getattr(my_module, 'CASEDATA')

    # network data
    nb = ds['xy'].shape[0]
    nl = ds['branch'].shape[0]
    nc = ds['cable'].shape[0]
    f, t, d, P1, Q1 = ds['branch'].T
    f = f.astype(int)
    t = t.astype(int)
    ibus = t
    Pd = np.zeros(nb)
    Qd = np.zeros(nb)
    Pd[ibus] = P1/1000
    Qd[ibus] = Q1/1000
    delta = ds['delta']
    area, r, x, Imax, ci = ds['cable'].T
    Imax = Imax/1000
    Smax = 3**0.5*ds['Vs']*Imax
    beta = 0.15*ds['alpha'] + 0.85*ds['alpha']**2
    Ci = ci*ds['g']
    xy = ds['xy']

    # incidence matrix
    Af = csr_matrix((np.ones(nl), (f, range(nl))), shape=(nb, nl))
    At = csr_matrix((np.ones(nl), (t, range(nl))), shape=(nb, nl))
    A = Af - At

    # load

    Pmin = Pd*(1 - delta/100)
    Pmax = Pd*(1 + delta/100)

    data = {
        'case_file': case_file,
        'Vs': ds['Vs'],
        'Vmin': ds['Vmin'],
        'f': f,
        't': t,
        'd': d,
        'A': A,
        'Ci': Ci,
        'nc': nc,
        'nl': nl,
        'nb': nb,
        'r': r,
        'x': x,
        'Pd': Pd,
        'Qd': Qd,
        'Pmin': Pmin,
        'Pmax': Pmax,
        'Smax': Smax,
        'beta': beta,
        'cl': ds['cl'],
        'delta': ds['delta'],
        'xy': xy,
    }

    return data


def deterministic(data):
    nb, nl = data['A'].shape
    nc = data['nc']
    # incidnce matrix without the first column representing slack bus
    A2 = data['A'][1:nb, :].todense()

    DP = pyo.ConcreteModel(name='Deterministc model')
    # SETS
    DP.lines = pyo.Set(initialize=range(nl))
    DP.buses = pyo.Set(initialize=range(nb))
    DP.loads = pyo.Set(initialize=range(1, nb))
    DP.cables = pyo.Set(initialize=range(nc))

    # PARAMETERS
    DP.P = pyo.Param(DP.lines, mutable=True)
    DP.Q = pyo.Param(DP.lines, mutable=True)
    DP.S = pyo.Param(DP.lines, mutable=True)
    # power flow in branch
    for i in DP.lines:
        DP.P[i] = -sum(np.linalg.inv(A2)[i, j] * data['Pd'][j+1]
                       for j in DP.lines)
        DP.Q[i] = -sum(np.linalg.inv(A2)[i, j] * data['Qd'][j+1]
                       for j in DP.lines)
        DP.S[i] = (DP.P[i]**2 + DP.Q[i]**2)**0.5
    # VARIABLES
    # cable status
    DP.v = pyo.Var(DP.lines, DP.cables, within=pyo.Binary)
    # square of bus voltages
    DP.W = pyo.Var(DP.buses, bounds=(data['Vmin']**2, None))
    # difference of square of voltages at both line ends
    DP.U = pyo.Var(DP.lines, bounds=(-data['Vs']**2, data['Vs']**2))

    # COSTS
    # line cost ($)
    DP.Cc = sum(sum(DP.v[i, j]*data['Ci'][j]
                for j in DP.cables)*data['d'][i] for i in DP.lines)
    # cost of energy losses ($)
    constant = 8760*data['cl']*data['beta']*1000/data['Vs']**2
    DP.Cl = sum(sum(DP.v[i, j]*data['r'][j] for j in DP.cables)
                * data['d'][i]*(DP.P[i]**2+DP.Q[i]**2)*constant for i in DP.lines)
    # Objective
    DP.obj = pyo.Objective(expr=DP.Cc + DP.Cl)

    # CONSTRAINTS
    # only one cable per line
    DP.oneCable = pyo.ConstraintList()
    for i in DP.lines:
        DP.oneCable.add(expr=sum(DP.v[i, j] for j in DP.cables) == 1)

    # supply bus voltage
    DP.supplyBus = pyo.Constraint(expr=DP.W[0] == data['Vs']**2)

    # line voltage equation
    DP.lineVoltage = pyo.ConstraintList()
    for i in DP.lines:
        DP.lineVoltage.add(expr=DP.U[i] == sum(
            data['A'][j, i]*DP.W[j] for j in DP.buses))

    # square voltage losses
    DP.voltageLosses = pyo.ConstraintList()
    for i in DP.lines:
        DP.voltageLosses.add(expr=DP.U[i] == 2*(DP.P[i]*sum(DP.v[i, j]*data['r'][j]
                             for j in DP.cables)+DP.Q[i]*sum(DP.v[i, j]*data['x'][j] for j in DP.cables))*data['d'][i])

    # max power flow in line
    DP.maxPowerFlow = pyo.ConstraintList()
    for i in DP.lines:
        DP.maxPowerFlow.add(expr=DP.S[i] <= sum(
            DP.v[i, j]*data['Smax'][j] for j in DP.cables))

    return DP


def subproblem(data, gama):
    nb, nl = data['A'].shape
    nc = data['nc']
    SP = pyo.ConcreteModel(name='Subproblem')

    # SETS
    SP.lines = pyo.Set(initialize=range(nl))
    SP.buses = pyo.Set(initialize=range(nb))
    SP.loads = pyo.Set(initialize=range(1, nb))
    SP.cables = pyo.Set(initialize=range(nc))

    # PARAMETERS
    SP.v = pyo.Param(SP.lines, SP.cables, mutable=True, default=1)

    # VARIABLES

    # active load demand
    SP.Pd = pyo.Var(SP.loads)
    # reactive power demands
    SP.Qd = pyo.Var(SP.loads)
    # linearize abs(Pd - Pd^{ref}) using new variable t
    SP.t = pyo.Var(SP.loads)
    # line active power flow
    # boundig P to the Smax of the biggest cable
    SP.P = pyo.Var(
        SP.lines, bounds=(0, data['Smax'][data['Smax'].shape[0]-1]))
    # line reactive power flow
    # boundig Q to the Smax of the biggest cable
    SP.Q = pyo.Var(
        SP.lines, bounds=(0, data['Smax'][data['Smax'].shape[0]-1]))
    # sqaure of bus voltages
    SP.W = pyo.Var(SP.buses, bounds=(data['Vmin']**2, None))
    # difference of square of voltages at both line ends
    SP.U = pyo.Var(SP.lines, bounds=(-data['Vs']**2, data['Vs']**2))

    # COSTS
    # cost of energy losses ($)
    constant = 8760*data['cl']*data['beta']*1000/data['Vs']**2
    SP.Cl = sum(sum(SP.v[i, j]*data['r'][j] for j in SP.cables)
                * data['d'][i]*(SP.P[i]**2+SP.Q[i]**2)*constant for i in SP.lines)
    # Objective
    SP.obj = pyo.Objective(expr=SP.Cl, sense=pyo.maximize)

    # CONSTRAINTS
    # supply bus voltage
    SP.supply = pyo.Constraint(expr=SP.W[0] == data['Vs']**2)

    # line voltage equation
    SP.lineVoltage = pyo.ConstraintList()
    for i in SP.lines:
        SP.lineVoltage.add(expr=SP.U[i] == sum(
            data['A'][j, i]*SP.W[j] for j in SP.buses))

    # square voltage losses U == 2*(P*(v*r*d) + Q*(v*x*d)) (v = const)
    SP.voltageLosses = pyo.ConstraintList()
    for i in SP.lines:
        SP.voltageLosses.add(
            expr=SP.U[i] == 2*(SP.P[i]*sum(SP.v[i, j]*data['r'][j] for j in SP.cables)*data['d'][i] + SP.Q[i]*sum(SP.v[i, j]*data['x'][j] for j in SP.cables)*data['d'][i]))

    # limits on line power flows
    SP.limitPQ = pyo.ConstraintList()
    # limits in active power flow
    for i in SP.lines:
        SP.limitPQ.add(SP.P[i] <= sum(SP.v[i, j]*data['Smax'][j]
                                      for j in SP.cables))
    # limits in reactive power flow
    for i in SP.lines:
        SP.limitPQ.add(SP.Q[i] <= sum(SP.v[i, j]*data['Smax'][j]
                                      for j in SP.cables))

    # limits on apparent power flow with piecewise linear terms to be used with CPLEX: P**2 + Q**2 <= Smax**2
    def parabola(SP, j, x):
        # we use j as the index for the constraint
        return x**2

    # piecewise constraints
    SP.limit_s_pwl = pyo.ConstraintList()
    # pwl variable for active and reactive power flow
    SP.P2 = pyo.Var(SP.lines)
    SP.Q2 = pyo.Var(SP.lines)
    # interpolation points
    pts = {}
    for i in SP.lines:
        s = sum(pyo.value(SP.v[i, j])*data['Smax'][j] for j in SP.cables)
        pts[i] = np.linspace(0, s, 10, endpoint=True).tolist()
    # constraints replacing P with pwl appriximation
    SP.pwl_p2 = pyo.Piecewise(SP.lines, SP.P2, SP.P,
                              pw_pts=pts, pw_constr_type='EQ', f_rule=parabola)
    # constraints replacing Q with pwl appriximation
    SP.pwl_q2 = pyo.Piecewise(SP.lines, SP.Q2, SP.Q,
                              pw_pts=pts, pw_constr_type='EQ', f_rule=parabola)
    # adding linearized constraints
    for i in SP.lines:
        SP.limit_s_pwl.add(
            SP.P2[i] + SP.Q2[i] <= sum(SP.v[i, j]*data['Smax'][j] for j in SP.cables)**2)

    # load balance for active power A * P == -Pd
    SP.loadBalance = pyo.ConstraintList()
    for i in SP.loads:
        SP.loadBalance.add(expr=sum(data['A'][i, j]*SP.P[j]
                                    for j in SP.lines) == -SP.Pd[i])

    # load balance for active power A * Q == -Qd
    for i in SP.loads:
        SP.loadBalance.add(
            expr=sum(data['A'][i, j]*SP.Q[j] for j in SP.lines) == -SP.Qd[i])

    # robust set
    SP.uncertaintySet = pyo.ConstraintList()
    for i in SP.loads:
        SP.uncertaintySet.add(SP.Pd[i] <= data['Pmax'][i])
        SP.uncertaintySet.add(SP.Pd[i] >= data['Pmin'][i])
        SP.uncertaintySet.add(SP.Qd[i] == data['Qd'][i]/data['Pd'][i]*SP.Pd[i])

    # linearize abs(Pd - Pd^{ref}) using new variable t
    # Pd^{ref} = data['Pd']
    # Pd^{Delta} = data['Pd'] - data['Pmin']

    SP.absPd = pyo.ConstraintList()
    for i in SP.loads:
        SP.absPd.add(SP.t[i] >= SP.Pd[i] - data['Pd'][i])
        SP.absPd.add(SP.t[i] >= data['Pd'][i] - SP.Pd[i])
    SP.gama = pyo.Constraint(expr=sum(
        SP.t[i] for i in SP.loads) <= gama*sum(data['Pd'][i]-data['Pmin'][i] for i in SP.loads))

    # LINE COST
    SP.Cc = sum(sum(SP.v[i, j]*data['Ci'][j]
                for j in SP.cables)*data['d'][i] for i in SP.lines)

    return SP


def master_problem(data):
    nb, nl = data['A'].shape
    nc = data['nc']
    MP = pyo.ConcreteModel(name='Master problem')

    # SETS
    MP.it = pyo.Set(initialize=[1])
    MP.lines = pyo.Set(initialize=range(nl))
    MP.buses = pyo.Set(initialize=range(nb))
    MP.loads = pyo.Set(initialize=range(1, nb))
    MP.cables = pyo.Set(initialize=range(nc))

    # PARAMETERS
    MP.Pd = pyo.Param(MP.loads, MP.it, mutable=True)
    MP.Qd = pyo.Param(MP.loads, MP.it, mutable=True)
    MP.P = pyo.Param(MP.lines, MP.it, mutable=True)
    MP.Q = pyo.Param(MP.lines, MP.it, mutable=True)
    MP.S = pyo.Param(MP.lines, MP.it, mutable=True)

    # VARIABLES
    # cable status
    MP.v = pyo.Var(MP.lines, MP.cables, within=pyo.Binary)
    # link with the subproblem
    MP.eta = pyo.Var(within=pyo.NonNegativeReals)
    # square of bus voltages
    MP.W = pyo.Var(MP.buses, MP.it, bounds=(data['Vmin']**2, None))
    # difference of square of voltages at both line ends
    MP.U = pyo.Var(MP.lines, MP.it, bounds=(-data['Vs']**2, data['Vs']**2))

    # COSTS
    # line cost ($)
    MP.Cc = sum(sum(MP.v[i, j]*data['Ci'][j]
                for j in MP.cables)*data['d'][i] for i in MP.lines)
    # Objective
    MP.obj = pyo.Objective(expr=MP.Cc + MP.eta)

    # CONSTRAINTS
    # only one cable per line
    MP.oneCable = pyo.ConstraintList()
    for i in MP.lines:
        MP.oneCable.add(expr=sum(MP.v[i, j] for j in MP.cables) == 1)
    # limit on eta
    MP.etaLimit = pyo.ConstraintList()
    # supply bus voltage
    MP.supplyBus = pyo.ConstraintList()
    # line voltage equation
    MP.lineVoltage = pyo.ConstraintList()
    # square voltage losses
    MP.voltageLosses = pyo.ConstraintList()
    # max power flow in line
    MP.maxPowerFlow = pyo.ConstraintList()

    return MP
