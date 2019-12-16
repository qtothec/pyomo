#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and 
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain 
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

from pyomo.common.config import (
    ConfigBlock, ConfigList, ConfigValue, In, NonNegativeFloat, NonNegativeInt,
    add_docstring_list, PositiveInt
)
from pyomo.common.fileutils import Executable
from pyomo.solver.base import OptSolver

class GUROBI(OptSolver):
    CONFIG = OptSolver.CONFIG()
    CONFIG.declare("solver_io", ConfigValue(
        default='lp',
        domain=In(['lp','nl','mps'])
    ))
    CONFIG.declare("executable", ConfigValue(
        default='gurobi.sh',
        domain=Executable,
    ))

    MAPPED_OPTIONS = OptSolver.MAPPED_OPTIONS()

    def solve(self, model, options=None, mapped_options=None, **config_options):
        """Solve a model""" + add_docstring_list("", GUROBI.CONFIG)

        options = self.OPTIONS(options)
        mapped_options = self.MAPPED_OPTIONS(mapped_options)
        config = self.CONFIG(config_options)

        
if __name__ == '__main__':
    a = GUROBI()
    print a.config.executable

    a.config.executable = 'emacs'
    print a.config.executable
