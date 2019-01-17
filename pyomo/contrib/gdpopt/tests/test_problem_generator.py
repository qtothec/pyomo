"""GDP problem generator"""
import random

import pyutilib.th as unittest

from pyomo.environ import ConcreteModel, RangeSet, Var, Integers, Binary, Objective, SolverFactory, \
    TransformationFactory
from pyomo.environ import TerminationCondition as tc
from pyomo.core.expr.current import exp, log
from pyomo.gdp import Disjunct


class TestProblemGenerator(unittest.TestCase):
    def test_generate_feasible_block(self):
        m = generate_feasible_block(ConcreteModel())
        build_block(m)

    def test_generate_model(self):
        m = generate_GDP_model()
        # SolverFactory('gdpopt').solve(m, mip_solver='gams', nlp_solver='gams')
        TransformationFactory('gdp.bigm').apply_to(m, bigM=500)
        SolverFactory('gams').solve(m, solver='baron', tee=True)


def generate_feasible_block(blk):
    max_iter = 25
    for retries in range(1, max_iter + 1):
        try:
            m = generate_block(ConcreteModel())
            build_block(m)
            m.objective = Objective(expr=m.x[1])
            # TODO suppress the warnings that I get from this.
            res = SolverFactory('gams').solve(m)
            res_tc = res.solver.termination_condition
            if res_tc in (tc.optimal, tc.feasible, tc.locallyOptimal):
                blk.var_data = m.var_data
                blk.bilinear_list = m.bilinear_list
                blk.exponential_list = m.exponential_list
                blk.log_list = m.log_list
                blk.poly_list = m.poly_list
                blk.frac_list = m.frac_list
                blk.linear_list = m.linear_list
                return blk
            else:
                # TODO use logger rather than print
                print("Termination condition of %s" % res_tc)
        except ValueError as e:
            # TODO use logger rather than print
            print("Feasible block generation attempt %s failed." % retries)
    raise RuntimeError("Unable to generatee feasible block in %s iterations." % max_iter)


def build_block(blk):
    n_vars, n_binary, n_integer = blk.var_data
    blk.vars = RangeSet(n_vars)
    blk.x = Var(blk.vars, bounds=(-100, 100), initialize=1)
    blk.binary_vars = RangeSet(n_binary)
    blk.integer_vars = RangeSet(n_integer)
    blk.y = Var(blk.binary_vars, domain=Binary)
    blk.z = Var(blk.integer_vars, domain=Integers)

    blk.bilinear = RangeSet(len(blk.bilinear_list))
    blk.exponential = RangeSet(len(blk.exponential_list))
    blk.logarithm = RangeSet(len(blk.log_list))
    blk.polynomial = RangeSet(len(blk.poly_list))
    blk.fractional = RangeSet(len(blk.frac_list))
    blk.linear = RangeSet(len(blk.linear_list))

    @blk.Constraint(blk.bilinear)
    def c_bilinear(blk, idx):
        A, i, j, k = blk.bilinear_list[idx - 1]
        return blk.x[i] <= A * blk.x[j] * blk.x[k]

    @blk.Constraint(blk.exponential)
    def c_exponential(blk, idx):
        A, B, i, j = blk.exponential_list[idx - 1]
        return blk.x[i] <= A * exp(B * blk.x[j])

    @blk.Constraint(blk.logarithm)
    def c_logarithm(blk, idx):
        A, B, i, j = blk.log_list[idx - 1]
        return blk.x[i] <= A * log(B * blk.x[j])

    @blk.Constraint(blk.polynomial)
    def c_polynomial(blk, idx):
        i, j = blk.poly_list[idx - 1][-2:]
        Ks = blk.poly_list[idx - 1][:-2]
        max_power = len(Ks) - 1
        return blk.x[i] <= sum(
            K * pow(10, max_power - p) * pow(blk.x[j], max_power - p)
            for p, K in enumerate(Ks))

    @blk.Constraint(blk.fractional)
    def c_fractional(blk, idx):
        A, B, i, j = blk.frac_list[idx - 1]
        return blk.x[i] <= A / (blk.x[j] + B)

    @blk.Constraint(blk.linear)
    def c_linear(blk, idx):
        A, B, i, j = blk.linear_list[idx - 1]
        return blk.x[i] <= A * blk.x[j] + B


def generate_GDP_model():
    n_disjunctions = 2
    disjuncts_per_disjunction = 2
    m = ConcreteModel()
    m.disjunctions = RangeSet(n_disjunctions)
    m.disjuncts = RangeSet(disjuncts_per_disjunction)
    m.disj = Disjunct(m.disjunctions, m.disjuncts)
    for disj in m.disj[...]:
        generate_feasible_block(disj)
        build_block(disj)

    # TODO improve naming convention
    @m.Disjunction(m.disjunctions)
    def c_disjunction(m, disjctn):
        return [m.disj[disjctn, disj] for disj in m.disjuncts]

    m.objective = Objective(expr=sum(m.disj[:, :].x[1]))

    m.data = {
        'n_disjunctions': n_disjunctions,
        'n_disjuncts': disjuncts_per_disjunction,
        'disjunctions': [
            [disj.data for disj in m.disj[disjctn, :]] for disjctn in m.disjunctions
        ]
    }
    return m


def generate_block(blk):
    # Variables
    n_vars = 10
    # TODO add support for binary and integer variables
    n_binary = 0
    n_integer = 0

    # Constraints
    n_bilinear = 3
    n_exponential = 1
    n_log = 1
    n_polynomial = 2
    poly_power = 3
    n_fractional = 1

    # Create model
    blk.var_data = (n_vars, n_binary, n_integer)

    """Bilinear constraints, of form:
    x[i] <= A x[j] x[k]
    A = random number
    
    m.bilinear_list = list(idx -> [A, i, j, k] for idx in m.bilinear)
    """
    blk.bilinear_list = _get_param_list(n_vars, n_bilinear, 1, 3)

    """Exponential constraints, of form:
    x[i] <= A * exp(B * x[j])
    A, B = random numbers
    
    m.exponential_list = list(idx -> [A, B, i, j] for idx in m.exponential)
    """
    blk.exponential_list = _get_param_list(n_vars, n_exponential, 2, 2)

    """Log function constraints, of form:
    x[i] <= A * log(B * x[j])
    A, B = random numbers
    
    m.log_list = list(idx -> [A, B, i, j] for idx in m.logarithm)
    """
    blk.log_list = _get_param_list(n_vars, n_log, 2, 2)

    """Polynomial constraints, of form:
    x[i] <= K0 * 10^-P * x[j]^P + K1 * 10^-(P-1) * x[j]^(P-1) + ... + K(P+1)
    P = maximum power
    K0, K1, ..., K(P+1) = random numbers
    
    m.poly_list = list(idx -> [K0, K1, ..., K(P+1), i, j] for idx in m.polynomial)
    """
    blk.poly_list = _get_param_list(n_vars, n_polynomial, poly_power + 1, 2)

    """Fractional constraints, of form:
    x[i] = A / (x[j] + B)
    A, B = random numbers
    
    m.frac_list = list(idx -> [A, B, i, j] for idx in m.fractional)
    """
    blk.frac_list = _get_param_list(n_vars, n_fractional, 2, 2)

    """Linear constraints, of form:
    x[i] <= A * x[j] + B"""
    # generate linear constraints so that everything is connected up to each other
    blk.used_var_idxs = set(
        [tup[-3:] for tup in blk.bilinear_list] +
        [tup[-2:] for tup in blk.exponential_list] +
        [tup[-2:] for tup in blk.log_list] +
        [tup[-2:] for tup in blk.poly_list] +
        [tup[-2:] for tup in blk.frac_list]
    )
    blk.unused_vars = set(range(1, n_vars + 1)) - blk.used_var_idxs
    vars_to_use = random.sample(blk.unused_vars, len(blk.unused_vars))
    blk.linear_list = []
    while len(vars_to_use) > 0:
        i, j = 0, 0
        num_vars_left = len(vars_to_use)
        if num_vars_left > 1:
            i, j = vars_to_use[:2]
            vars_to_use = vars_to_use[2:]
        elif num_vars_left == 1:
            i = vars_to_use[0]
            j = random.sample(blk.used_var_idxs, 1)
            vars_to_use = []
        blk.linear_list.append((
            random.uniform(-10, 10), random.uniform(-10, 10), i, j))

    return blk


def _get_param_list(tot_vars, n, n_rnd, n_vars):
    return list(
        tuple(random.uniform(-10, 10) for _ in range(n_rnd)) + tuple(_get_n_var_idx(tot_vars, n_vars))
        for _ in range(n))


def _get_n_var_idx(tot_vars, n):
    return random.sample(range(1, tot_vars + 1), n)


if __name__ == "__main__":
    unittest.main()
