#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and 
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain 
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________


class GUROBI(OptSolver):
    CONFIG = OptSolver.CONFIG()
    #CONFIG.declare(...)

    MAPPED_OPTIONS = OptSolver.MAPPED_OPTIONS()

    def solve(self, model, options=None, mapped_options=None, **config_options):
        """Solve a model""" + add_docstring_list("", GUROBI.CONFIG)

        options = self.OPTIONS(options)
        mapped_options = self.MAPPED_OPTIONS(mapped_options)
        config = self.CONFIG(config_options)

        
