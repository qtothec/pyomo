#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and 
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain 
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________


from pyomo.common.config import ConfigBlock, ConfigValue, NonNegativeFloat
from pyomo.common.errors import DeveloperError
from pyomo.common.deprecation import deprecated

from pyomo.opt.results import SolverResults

class Solver(object):
    """A generic optimization solver"""

    CONFIG = ConfigBlock()
    CONFIG.declare('timelimit', ConfigValue(
        default=None,
        domain=NonNegativeFloat,
    ))
    CONFIG.declare('keepfiles', ConfigValue(
        default=False,
        domain=bool,
    ))
    CONFIG.declare('tee', ConfigValue(
        default=False,
        domain=bool,
    ))
    CONFIG.declare('load_solution', ConfigValue(
        default=True,
        domain=bool,
    ))


    def __init__(self, **kwds):
        self.config = self.CONFIG()
        self.options = ConfigBlock(implicit=True)

    def available(self):
        raise DeveloperError(
            "Derived Solver class %s failed to implement available()"
            % (self.__class__.__name__,))

    def license_status(self):
        raise DeveloperError(
            "Derived Solver class %s failed to implement license_status()"
            % (self.__class__.__name__,))
        
    def version(self):
        """
        Returns a tuple describing the solver version.
        """
        raise DeveloperError(
            "Derived Solver class %s failed to implement version()"
            % (self.__class__.__name__,))

    def solve(self, model, options=None, **config_options):
        raise DeveloperError(
            "Derived Solver class %s failed to implement solve()"
            % (self.__class__.__name__,))

    @deprecated("Casting a solver to bool() is deprecated.  Use available()",
                version='TBD')
    def __bool__(self):
        return self.available()


class MIPSolver(Solver):
    CONFIG = Solver.CONFIG()
    CONFIG.declare('mipgap', ConfigValue(
        default=None,
        domain=NonNegativeFloat,
    ))
    CONFIG.declare('relax_integrality', ConfigValue(
        default=False,
        domain=bool,
    ))
    
