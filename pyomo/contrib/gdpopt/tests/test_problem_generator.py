"""GDP problem generator"""
import json
import os
import random

import pyutilib.th as unittest
from pyutilib.misc import Container

from pyomo.environ import ConcreteModel, RangeSet, Var, Integers, Binary, Objective, SolverFactory, \
    TransformationFactory, value
from pyomo.environ import TerminationCondition as tc
from pyomo.core.expr.current import exp, log


class TestProblemGenerator(unittest.TestCase):
    def test_generate_model(self):
        m = build_GDP_model(generate_GDP_model_data())
        TransformationFactory('gdp.bigm').apply_to(m, bigM=500)
        SolverFactory('gams').solve(m, solver='baron', tee=True)

    def test_read_write_model(self):
        blk_data = generate_GDP_model_data()
        with open("test.json", 'w') as output_file:
            json.dump(blk_data, output_file)
        with open("test.json") as input_file:
            recover = json.load(input_file)
        recovered_blk_data = Container()
        recovered_blk_data.update(recover)
        m = build_GDP_model(recovered_blk_data)
        TransformationFactory('gdp.bigm').apply_to(m, bigM=500)
        SolverFactory('gams').solve(m, solver='baron', tee=True)
        os.remove("test.json")


def next_feasible_block_data():
    """Generates initial feasible block data."""
    max_iter = 25
    for retries in range(1, max_iter + 1):
        print('Attempt {0}'.format(retries))
        try:
            blk_data = generate_basic_block_data()
            m = build_basic_block(ConcreteModel(), blk_data)
            m.objective = Objective(expr=m.x[0])
            # TODO suppress the warnings that I get from this.
            res = SolverFactory('gams').solve(m)
            res_tc = res.solver.termination_condition
            if res_tc in (tc.optimal, tc.feasible, tc.locallyOptimal):
                # need to do additional check if model objective value is at bound
                if value(m.x[0]) == m.x[0].lb:
                    # block is unbounded. Regenerate.
                    print("Unbounded.")
                    continue
                return blk_data
            else:
                # TODO use logger rather than print
                print("Termination condition of %s" % res_tc)
        except ValueError as e:
            # TODO use logger rather than print
            print("Feasible block generation attempt %s failed." % retries)
    raise RuntimeError("Unable to generate feasible block in %s iterations." % max_iter)


class ModelBlockData(Container):
    """Data class to hold model block information

    Attributes:
        - version: data class version
        - variables: dictionary of variable information
            - continuous: list of tuples with variable information: (lb, ub, init_val)
            - binary: list of tuples with variable information: (lb, ub, init_val)
            - integer: list of tuples with variable information: (lb, ub, init_val)
        - constraints: dictionary of constraint information
            - bilinear: list of tuples (A, i, j, k) --> x[i] <= A x[j] x[k]
            - exponential: list of tuples (A, B, i, j) --> x[i] <= A exp(B x[j])
            - logarithm: list of tuples (A, B, i, j) --> x[i] <= log(B x[j])
            - polynomial: list of tuples (K0, K1, ..., K(P+1), i, j)
                - P = maximum power
                - x[i] <= K0 10^-P + x[j]^P + K1 10^(-(P-1)) x[j]^(P-1) + ... + K[P+1]
            - fractional: list of tuples (A, B, i, j) --> x[i] <= A / (x[j] + B)
            - linear: (A, B, i, j) --> x[i] <= A x[j] + B
        - disjuncts: list of ModelBlockData objects
        - disjunctions: list of tuples (D0, D1, ..., DN) --> D0 V D1 V ... V DN
            - D0, D1, ..., DN are disjunct indices
    """
    version = 1


def generate_GDP_model_data():
    n_disjunctions = 2
    max_disjuncts_per_disjunction = 4

    blk_data = next_feasible_block_data()
    num_disjuncts = [random.randint(2, max_disjuncts_per_disjunction) for _ in range(n_disjunctions)]
    counter = 0
    for i, num in enumerate(num_disjuncts):
        num_disjuncts[i] = tuple(range(counter, counter + num))
        counter += num
    blk_data.disjunctions = num_disjuncts
    blk_data.disjuncts = [next_feasible_block_data() for _ in range(counter)]

    return blk_data


def generate_basic_block_data():
    """Generate the initial block data."""
    n_continuous = 10
    n_positive = 3

    # constraints
    n_bilinear = 3
    n_exponential = 1
    n_log = 1
    n_polynomial = 2
    poly_power = 3
    n_fractional = 1

    used_vars = set()

    def sample_vars(n):
        sample = tuple(random.sample(range(n_continuous), n))
        used_vars.update(sample)
        return sample

    positive_vars = sample_vars(n_positive)

    def sample_positive_vars(n):
        sample = tuple(random.sample(positive_vars, n))
        used_vars.update(sample)
        return sample

    model_data = ModelBlockData()
    model_data.variables = Container(
        continuous=[(-100, 100, 1) for _ in range(n_continuous)],
        binary=[],
        integer=[],
    )
    for v_idx in positive_vars:
        model_data.variables.continuous[v_idx] = (0.001, 100, 1)
    model_data.constraints = Container(
        bilinear=[(random.uniform(-10, 10),) + sample_vars(3)
                  for _ in range(n_bilinear)],
        exponential=[(random.uniform(-10, 10), random.uniform(-2, 2),) + sample_vars(2)
                     for _ in range(n_exponential)],
        logarithm=[(random.uniform(-10, 10), random.uniform(0.001, 2),) + sample_positive_vars(2)
                   for _ in range(n_log)],
        polynomial=[tuple(random.uniform(-10, 10) for _ in range(poly_power + 1)) + sample_vars(2)
                    for _ in range(n_polynomial)],
        fractional=[(random.uniform(-10, 10), random.uniform(-10, 10),) + sample_vars(2)
                    for _ in range(n_fractional)],
        linear=[],  # defined below
    )

    unused_vars = set(range(n_continuous)) - used_vars
    vars_to_use = random.sample(unused_vars, len(unused_vars))
    num_vars_left = len(vars_to_use)
    while num_vars_left > 0:
        i, j = 0, 0
        if num_vars_left > 1:
            i, j = vars_to_use[:2]
            vars_to_use = vars_to_use[2:]
        elif num_vars_left == 1:
            i = vars_to_use[0]
            j = random.choice(tuple(used_vars))
            vars_to_use = []
        model_data.constraints['linear'].append((random.uniform(-10, 10), random.uniform(-10, 10), i, j))
        num_vars_left = len(vars_to_use)

    return model_data


def build_basic_block(blk, blk_data):
    blk._prob_data = blk_data
    blk.continuous_vars = RangeSet(0, len(blk_data.variables['continuous']) - 1)
    blk.x = Var(blk.continuous_vars)
    for i, var_data in enumerate(blk_data.variables['continuous']):
        blk.x[i].setlb(var_data[0])
        blk.x[i].setub(var_data[1])
        blk.x[i].set_value(var_data[2])
    blk.binary_vars = RangeSet(0, len(blk_data.variables['binary']) - 1)
    blk.y = Var(blk.binary_vars, domain=Binary)
    for i, var_data in enumerate(blk_data.variables['binary']):
        blk.y[i].setlb(var_data[0])
        blk.y[i].setub(var_data[1])
        blk.y[i].set_value(var_data[2])
    blk.integer_vars = RangeSet(0, len(blk_data.variables['integer']) - 1)
    blk.z = Var(blk.integer_vars, domain=Integers)
    for i, var_data in enumerate(blk_data.variables['integer']):
        blk.z[i].setlb(var_data[0])
        blk.z[i].setub(var_data[1])
        blk.z[i].set_value(var_data[2])

    blk.bilinear = RangeSet(0, len(blk_data.constraints['bilinear']) - 1)
    blk.exponential = RangeSet(0, len(blk_data.constraints['exponential']) - 1)
    blk.logarithm = RangeSet(0, len(blk_data.constraints['logarithm']) - 1)
    blk.polynomial = RangeSet(0, len(blk_data.constraints['polynomial']) - 1)
    blk.fractional = RangeSet(0, len(blk_data.constraints['fractional']) - 1)
    blk.linear = RangeSet(0, len(blk_data.constraints['linear']) - 1)

    @blk.Constraint(blk.bilinear)
    def c_bilinear(blk, idx):
        A, i, j, k = blk_data.constraints['bilinear'][idx]
        return blk.x[i] <= A * blk.x[j] * blk.x[k]

    @blk.Constraint(blk.exponential)
    def c_exponential(blk, idx):
        A, B, i, j = blk_data.constraints['exponential'][idx]
        return blk.x[i] <= A * exp(B * blk.x[j])

    @blk.Constraint(blk.logarithm)
    def c_logarithm(blk, idx):
        A, B, i, j = blk_data.constraints['logarithm'][idx]
        return blk.x[i] <= A * log(B * blk.x[j])

    @blk.Constraint(blk.polynomial)
    def c_polynomial(blk, idx):
        i, j = blk_data.constraints['polynomial'][idx][-2:]
        Ks = blk_data.constraints['polynomial'][idx][:-2]
        max_power = len(Ks) - 1
        return blk.x[i] <= sum(
            K * pow(10, max_power - p) * pow(blk.x[j], max_power - p)
            for p, K in enumerate(Ks))

    @blk.Constraint(blk.fractional)
    def c_fractional(blk, idx):
        A, B, i, j = blk_data.constraints['fractional'][idx]
        return blk.x[i] <= A / (blk.x[j] + B)

    @blk.Constraint(blk.linear)
    def c_linear(blk, idx):
        A, B, i, j = blk_data.constraints['linear'][idx]
        return blk.x[i] <= A * blk.x[j] + B

    return blk


def build_GDP_model(model_data):
    m = ConcreteModel()
    m.disjunctions = RangeSet(0, len(model_data.disjunctions) - 1)
    m.disjuncts = RangeSet(0, len(model_data.disjuncts) - 1)

    @m.Disjunct(m.disjuncts)
    def d(disj, idx):
        return build_basic_block(disj, model_data.disjuncts[idx])

    @m.Disjunction(m.disjunctions)
    def c_disjunction(m, disjctn):
        return [m.d[disjunct_idx] for disjunct_idx in model_data.disjunctions[disjctn]]

    # Objective is currently to minimize the first variable in each disjunct.
    # TODO come up with a more intelligent objective function to link the different disjuncts
    m.objective = Objective(expr=sum(m.d[:].x[0]))

    return m


if __name__ == "__main__":
    unittest.main()
