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
from pyomo.solver.base import OptSolver

class GurobiSolver(OptSolver):
    CONFIG = OptSolver.CONFIG()
    CONFIG.declare("solver_io", ConfigValue(
        default='lp',
        domain=In(['lp','nl','mps'])
    ))

    def __new__(cls, **kwds):
        config = self.CONFIG()
        config.solver_io = kwds.pop('solver_io', config.solver_io)
        if config.solver_io == 'lp':
            return GurobiSolver_LP(*args, **kwds)
        if config.solver_io == 'nl':
            return GurobiSolver_NL(*args, **kwds)
        if config.solver_io == 'mps':
            return GurobiSolver_MPS(*args, **kwds)
        # For the direct / persistent solvers, they are implemented in
        # other modules.  To simplify imports, we will defer to the
        # SolverFactory
        if config.solver_io == 'persistent':
            return SolverFactory('gurobi_persistent', **kwds)
        if config.solver_io in ('direct', 'python'):
            return SolverFactory('gurobi_direct', **kwds)

class GurobiSolver_LP(GurobiSolver):
    CONFIG = GurobiSolver.CONFIG()
    CONFIG.declare("executable", ConfigValue(
        default='gurobi.bat' if sys.platform == 'win32' else 'gurobi.sh' ,
        domain=Executable,
    ))

    MAPPED_OPTIONS = OptSolver.MAPPED_OPTIONS()

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
        """Solve a model""" + add_docstring_list("", GUROBI.CONFIG)

        options = self.OPTIONS(options)
        mapped_options = self.MAPPED_OPTIONS(mapped_options)
        config = self.CONFIG(config_options)

        try:
            return self._apply_solver(model, options, mapped_options, config)
        finally:
            self.cleanup(model, options, mapped_options, config)


    def _apply_solver(self, model, options, mapped_options, config):
        if not config.problemfile is None:
            config.problemfile = TempfileManager.create_tempfile(
                suffix='.pyomo.%s' % (config.solver_io,))
        if not config.logfile:
            config.logfile = TempfileManager.create_tempfile(
                suffix='.gurobi.log')
        if not config.solnfile:
            config.solnfile = TempfileManager.create_tempfile(
                suffix='.gurobi.txt')
        
        writer = WriterFactory(config.solver_io)
        fname, symbol_map = writer(
            model=model,
            output_filename=config.problemfile.path(),
            solver_capability=lambda x: True,
            io_options=config.io_options
        )
        assert fname == str(config.problemfile)

        
               



import logging
logging.getLogger('pyomo.common').setLevel(logging.INFO)

if __name__ == '__main__':
    a = GUROBI()
    print a.config.executable

    a.config.executable = 'emacs'
    print a.config.executable
