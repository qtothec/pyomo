#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

import os
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle

from pyutilib.common import ApplicationError, WindowsError
from pyutilib.services import TempfileManager
from pyutilib.subprocess import run

from pyomo.common.config import (
    ConfigBlock, ConfigList, ConfigValue, add_docstring_list, In, Path,
)
from pyomo.common.fileutils import Executable, this_file_dir
from pyomo.solver.base import MIPSolver, SolverResults
from pyomo.writer.cpxlp import ProblemWriter_cpxlp

class GurobiSolver(MIPSolver):
    CONFIG = MIPSolver.CONFIG()

    def __new__(cls, **kwds):
        if cls != GurobiSolver:
            return super(GurobiSolver, cls).__new__(cls, **kwds)

        solver_io = kwds.pop('solver_io', 'lp')
        if solver_io == 'lp':
            return GurobiSolver_LP(**kwds)
        elif solver_io == 'nl':
            return GurobiSolver_NL(**kwds)
        elif solver_io == 'mps':
            return GurobiSolver_MPS(**kwds)
        # For the direct / persistent solvers, they are implemented in
        # other modules.  To simplify imports, we will defer to the
        # SolverFactory
        elif solver_io == 'persistent':
            return SolverFactory('gurobi_persistent', **kwds)
        elif solver_io in ('direct', 'python'):
            return SolverFactory('gurobi_direct', **kwds)
        else:
            raise ValueError("Invalid solver_io for GurobiSolver: %s"
                             % (solver_io,))


class GurobiSolver_LP(GurobiSolver):
    CONFIG = GurobiSolver.CONFIG()
    CONFIG.declare("executable", ConfigValue(
        default='gurobi.bat' if sys.platform == 'win32' else 'gurobi.sh' ,
        domain=Executable,
    ))
    CONFIG.declare("problemfile", ConfigValue(
        default=None,
        domain=Path(),
    ))
    CONFIG.declare("logfile", ConfigValue(
        default=None,
        domain=Path(),
    ))
    CONFIG.declare("solnfile", ConfigValue(
        default=None,
        domain=Path(),
    ))
    CONFIG.declare("warmstart_file", ConfigValue(
        default=None,
        domain=Path(),
    ))

    CONFIG.inherit_from(ProblemWriter_cpxlp.CONFIG, skip={
        'allow_quadratic_objective', 'allow_quadratic_constraints',
        'allow_sos1', 'allow_sos2'})

    def available(self):
        return self.config.executable.path() is not None

    def license_status(self):
        """
        Runs a check for a valid Gurobi license using the
        given executable (default is 'gurobi_cl'). All
        output is hidden. If the test fails for any reason
        (including the executable being invalid), then this
        function will return False.
        """
        if not self.available():
            return False
        gurobi_cl = os.path.join(os.path.dirname(self.config.executable),
                                 'gurobi_cl')
        try:
            rc = subprocess.call([gurobi_cl, "--license"],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
        except OSError:
            rc = 1
        return rc == 0


    def solve(self, model, options=None, **config_options):
        """Solve a model""" + add_docstring_list("", GurobiSolver_LP.CONFIG)

        options = self.options(options)
        config = self.config(config_options)

        try:
            TempfileManager.push()
            return self._apply_solver(model, options, config)
        finally:
            # finally, clean any temporary files registered with the
            # temp file manager, created populated *directly* by this
            # plugin.
            TempfileManager.pop(remove=not config.keepfiles)


    def _apply_solver(self, model, options, config):
        if not config.problemfile:
            config.problemfile = TempfileManager.create_tempfile(
                suffix='.pyomo.lp')
        if not config.logfile:
            config.logfile = TempfileManager.create_tempfile(
                suffix='.gurobi.log')
        if not config.solnfile:
            config.solnfile = TempfileManager.create_tempfile(
                suffix='.gurobi.txt')
        
        # Gurobi can support certain problem features
        writer_config = ProblemWriter_cpxlp.CONFIG()
        writer_config.allow_quadratic_objective = True
        writer_config.allow_quadratic_constraints = True
        writer_config.allow_sos1 = True
        writer_config.allow_sos2 = True
        # Copy over the relevant values from the solver config
        # (skip_implicit alloes the config to have additional fields
        # that are ignored)
        writer_config.set_value(config, skip_implicit=True)
        fname, symbol_map = ProblemWriter_cpxlp()(
            model=model,
            output_filename=config.problemfile,
            io_options=writer_config
        )
        assert fname == str(config.problemfile)

        # Handle mapped options
        mipgap = config.mipgap
        if mipgap is not None:
            options['MIPGap'] = mipgap
        options['LogFile'] = config.logfile

        # Extract the suffixes
        suffixes = []

        # Run Gurobi
        data = pickle.dumps(
            ( config.problemfile,
              config.solnfile,
              { 'warmstart_file': config.warmstart_file,
                'relax_integrality': config.relax_integrality, },
              options.value(),
              suffixes))
        timelim = config.timelimit
        if timelim:
            timelim + min(max(1, 0.01*self._timelim), 100)
        cmd = [ str(config.executable),
                os.path.join(this_file_dir(), 'GUROBI_RUN.py') ]
        try:
            rc, log = run(cmd, stdin=data, timelimit=timelim, tee=config.tee)
        except WindowsError:
            raise ApplicationError(
                'Could not execute the command: %s\tError message: %s'
                % (cmd, sys.exc_info()[1]))
        sys.stdout.flush()

        # Read in the results
        result_data = None
        with open(config.solnfile, 'rb') as SOLN:
            try:
                result_data = pickle.load(SOLN)
            except ValueError:
                logger.error(
                    "GurobiSolver_LP: no results data returned from the "
                    "Gurobi subprocess.  Look at the solver log for more "
                    "details (re-run with 'tee=True' to see the solver log.")
                raise

        results = SolverResults()
        results.problem.update(result_data['problem'])
        results.solver.update(result_data['solver'])
        results.solver.name = 'gurobi_lp'

        if not config.load_solution:
            raise RuntimeError("TODO")
        elif result_data['solution']['points']:
            _sol = result_data['solution']['points'][0]
            X = _sol['X']
            for i, vname in enumerate(_sol['VarName']):
                v = symbol_map.getObject(vname)
                v.value = X[i]

        return results
        

if __name__ == '__main__':
    solver = GurobiSolver_LP()
    from pyomo.environ import *
    m = ConcreteModel()
    m.x = Var(bounds=(0,10))
    m.y = Var(bounds=(0,5))
    m.con1 = Constraint(expr=2*m.x+m.y <= 5)
    m.con2 = Constraint(expr=m.x+2*m.y <= 6)
    m.obj = Objective(expr=m.x+m.y, sense=maximize)

    print solver.solve(m, tee=True)
    m.pprint()

    
