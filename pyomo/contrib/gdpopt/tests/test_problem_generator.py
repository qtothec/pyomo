"""GDP problem generator"""
import random

import pyutilib.th as unittest

from pyomo.core import ConcreteModel, RangeSet, Var, Integers, Binary, Objective
from pyomo.core.expr.current import exp, log


class TestProblemGenerator(unittest.TestCase):
    def test_build_model(self):
        m = generate_block()
        build_model(m)
        m.pprint()


def generate_model():
    n_disjunctions = 2
    disjuncts_per_disjunction = 2


def generate_block():
    # Variables
    n_vars = 10
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
    m = ConcreteModel()
    m.vars = RangeSet(n_vars)
    m.x = Var(m.vars, bounds=(-100, 100))
    m.binary_vars = RangeSet(n_binary)
    m.integer_vars = RangeSet(n_integer)
    m.y = Var(m.binary_vars, domain=Binary)
    m.z = Var(m.integer_vars, domain=Integers)
    m.objective = Objective(expr=m.x[1])

    """Bilinear constraints, of form:
    x[i] <= A x[j] x[k]
    A = random number
    
    m.bilinear_list = list(idx -> [A, i, j, k] for idx in m.bilinear)
    """
    m.bilinear = RangeSet(n_bilinear)
    m.bilinear_list = _get_param_list(m, len(m.bilinear), 1, 3)

    """Exponential constraints, of form:
    x[i] <= A * exp(B * x[j])
    A, B = random numbers
    
    m.exponential_list = list(idx -> [A, B, i, j] for idx in m.exponential)
    """
    m.exponential = RangeSet(n_exponential)
    m.exponential_list = _get_param_list(m, len(m.exponential), 2, 2)

    """Log function constraints, of form:
    x[i] <= A * log(B * x[j])
    A, B = random numbers
    
    m.log_list = list(idx -> [A, B, i, j] for idx in m.logarithm)
    """
    m.logarithm = RangeSet(n_log)
    m.log_list = _get_param_list(m, len(m.logarithm), 2, 2)

    """Polynomial constraints, of form:
    x[i] <= K0 * 10^-P * x[j]^P + K1 * 10^-(P-1) * x[j]^(P-1) + ... + K(P+1)
    P = maximum power
    K0, K1, ..., K(P+1) = random numbers
    
    m.poly_list = list(idx -> [K0, K1, ..., K(P+1), i, j] for idx in m.polynomial)
    """
    m.polynomial = RangeSet(n_polynomial)
    m.poly_list = _get_param_list(m, len(m.polynomial), poly_power + 1, 2)

    """Fractional constraints, of form:
    x[i] = A / (x[j] + B)
    A, B = random numbers
    
    m.frac_list = list(idx -> [A, B, i, j] for idx in m.fractional)
    """
    m.fractional = RangeSet(n_fractional)
    m.frac_list = _get_param_list(m, len(m.fractional), 2, 2)

    """Linear constraints, of form:
    x[i] <= A * x[j] + B"""
    # generate linear constraints so that everything is connected up to each other
    m.used_var_idxs = set(
        [tup[-3:] for tup in m.bilinear_list] +
        [tup[-2:] for tup in m.exponential_list] +
        [tup[-2:] for tup in m.log_list] +
        [tup[-2:] for tup in m.poly_list] +
        [tup[-2:] for tup in m.frac_list]
    )
    m.unused_vars = set(m.vars - m.used_var_idxs)
    vars_to_use = random.sample(m.unused_vars, len(m.unused_vars))
    m.linear_list = []
    while len(vars_to_use) > 0:
        i, j = 0, 0
        num_vars_left = len(vars_to_use)
        if num_vars_left > 1:
            i, j = vars_to_use[:2]
            vars_to_use = vars_to_use[2:]
        elif num_vars_left == 1:
            i = vars_to_use[0]
            j = random.sample(m.used_var_idxs, 1)
            vars_to_use = []
        m.linear_list.append((
            random.uniform(-10, 10), random.uniform(-10, 10), i, j))
    m.linear = RangeSet(len(m.linear_list))

    return m


def build_model(m):
    @m.Constraint(m.bilinear)
    def c_bilinear(m, idx):
        A, i, j, k = m.bilinear_list[idx - 1]
        return m.x[i] <= A * m.x[j] * m.x[k]

    @m.Constraint(m.exponential)
    def c_exponential(m, idx):
        A, B, i, j = m.exponential_list[idx - 1]
        return m.x[i] <= A * exp(B * m.x[j])

    @m.Constraint(m.logarithm)
    def c_logarithm(m, idx):
        A, B, i, j = m.log_list[idx - 1]
        return m.x[i] <= A * log(B * m.x[j])

    @m.Constraint(m.polynomial)
    def c_polynomial(m, idx):
        i, j = m.poly_list[idx - 1][-2:]
        Ks = m.poly_list[idx - 1][:-2]
        max_power = len(Ks) - 1
        return m.x[i] <= sum(
            K * pow(10, max_power - p) * pow(m.x[j], max_power - p)
            for p, K in enumerate(Ks))

    @m.Constraint(m.fractional)
    def c_fractional(m, idx):
        A, B, i, j = m.frac_list[idx - 1]
        return m.x[i] <= A / (m.x[j] + B)

    @m.Constraint(m.linear)
    def c_linear(m, idx):
        A, B, i, j = m.linear_list[idx - 1]
        return m.x[i] <= A * m.x[j] + B


def _get_param_list(m, n, n_rnd, n_vars):
    return list(
        tuple(random.uniform(-10, 10) for _ in range(n_rnd)) + tuple(_get_n_var_idx(m, n_vars))
        for _ in range(n))

def _get_n_var_idx(m, n):
    return random.sample(range(1, len(m.vars) + 1), n)


"""Nonlinear inequalities
One term, two term
Bilinear
Exponential
Log function
Polynomial (varying order)
Fractional terms"""


if __name__ == "__main__":
    unittest.main()
