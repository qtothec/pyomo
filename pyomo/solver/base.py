#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and 
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain 
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________


from pyomo.common.config import ConfigBlock
from pyomo.common.errors import DeveloperError
from pyomo.common.deprecation import deprecated

class OptSolver(object):
    """A generic optimization solver"""

    CONFIG = ConfigBlock()

    MAPPED_OPTIONS = ConfigBlock()

    def __init__(self, **kwds):
        self.config = OptSolver.CONFIG()
        self.mapped_options = OptSolver.MAPPED_OPTIONS()
        self.options = ConfigBlock(implicit=True)

    def available(self):
        raise DeveloperError(
            "Derived OptSolver class %s failed to implement available()"
            % (self.__class__.__name__,))

    def license_status(self):
        raise DeveloperError(
            "Derived OptSolver class %s failed to implement license_status()"
            % (self.__class__.__name__,))
        
    def version(self):
        """
        Returns a tuple describing the solver version.
        """
        raise DeveloperError(
            "Derived OptSolver class %s failed to implement version()"
            % (self.__class__.__name__,))

    def solve(self, model, options=None, mapped_options=None, **config_options):
        raise DeveloperError(
            "Derived OptSolver class %s failed to implement solve()"
            % (self.__class__.__name__,))

    @deprecated("Casting a solver to bool() is deprecated.  Use available()",
                version='TBD')
    def __bool__(self):
        return self.available()
