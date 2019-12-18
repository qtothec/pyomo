#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

import sys

from pyutilib.services import TempfileManager

from pyomo.common.config import (
    ConfigBlock, ConfigList, ConfigValue, add_docstring_list, In,
)
from pyomo.common.fileutils import Executable
from pyomo.solver.base import Solver
from pyomo.writer.cpxlp import ProblemWriter_cpxlp

class GurobiSolver(Solver):
    CONFIG = Solver.CONFIG()

    MAPPED_OPTIONS = Solver.MAPPED_OPTIONS()

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
    CONFIG.inherit_from(ProblemWriter_cpxlp.CONFIG, skip={
        'allow_quadratic_objective', 'allow_quadratic_constraints',
        'allow_sos1', 'allow_sos2'})

    MAPPED_OPTIONS = GurobiSolver.MAPPED_OPTIONS()

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


    def solve(self, model, options=None, mapped_options=None, **config_options):
        """Solve a model""" + add_docstring_list("", GurobiSolver_LP.CONFIG)

        options = self.options(options)
        mapped_options = self.mapped_options(mapped_options)
        config = self.config(config_options)

        try:
            return self._apply_solver(model, options, mapped_options, config)
        finally:
            self.cleanup(model, options, mapped_options, config)


    def _apply_solver(self, model, options, mapped_options, config):
        if not config.problemfile is None:
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
        writer_config.allow_quadratic_objective.set_default(True)
        writer_config.allow_quadratic_constrinat.set_default(True)
        writer_config.allow_sos1.set_default(True)
        writer_config.allow_sos2.set_default(True)
        # Copy over the relevant values from the solver config
        # (skip_implicit alloes the config to have additional fields
        # that are ignored)
        writer.set_value(config, skip_implicit=True)
        fname, symbol_map = ProblemWriter_cpxlp(
            model=model,
            output_filename=config.problemfile.path(),
            io_options=writer_config
        )
        assert fname == str(config.problemfile)

        # Handle mapped options
        mipgap = mapped_options.mipgap
        if mipgap is not None:
            options['MIPGap'] = mipgap

        # Extract the suffixes
        

        # Run Gurobi
        data = json.dumps(
            ( config.problemfile,
              config.solnfile,
              { 'warmstart_file': config.warmstart_file,
                'relax_integrality': mapped_options.relax_integrality, },
              options,
              suffixes))
        timelim = config.timelimit
        if timelim:
            timelim + min(max(1, 0.01*self._timelim), 100)
        cmd = [ config.executable,
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
        with open(config.solnfile, 'r') as SOLN:
            try:
                result_data = json.load(SOLN)
            except ValueError:
                logger.error(
                    "GurobiSolver_LP: no results data returned from the
                    Gurobi subprocess.  Look at the solver log for more
                    details (re-run with 'tee=True' to see the solver log.")
                raise

        
            

if __name__ == '__main__':
    a = GurobiSolver_LP()
    a.CONFIG.display()

    
