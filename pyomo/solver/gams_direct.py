import logging

from pyutilib.misc.config import ConfigValue
from pyutilib.misc.timing import TicTocTimer

from pyomo.common.config import add_docstring_list
from pyomo.opt import SolverFactory
from pyomo.solver.base import Solver

logger = logging.getLogger('pyomo.solvers')


@SolverFactory.register('new_gams_direct', doc="Direct python interface to GAMS")
class GAMSDirect(Solver):
    """Direct interface to GAMS"""
    CONFIG = Solver.CONFIG()

    CONFIG.declare('symbolic_solver_labels', ConfigValue(
        default=False, domain=bool,
        doc="If True, the GAMS variable and constraint names will "
        "closely match those of the pyomo variables and constraints."
    ))
    CONFIG.declare('stream_solver', ConfigValue(
        default=False, domain=bool,
        doc="If True, show the GAMS output"
    ))

    __doc__ = add_docstring_list(__doc__, CONFIG)

    def __init__(self):
        super(GAMSDirect, self).__init__()

        try:
            import gams
            self._gams_module = gams
            self._python_api_exists = True
        except Exception as e:
            logger.warning("Import of GAMS Python Interface failed - GAMS message=" + str(e) + "\n")
            self._python_api_exists = False

    def available(self):
        return self._python_api_exists

    def license_status(self):
        if not self.available():
            return False
        # TODO figure out how to check license status
        return True

    def solve(self, model, options=None, **config_options):
        """Solve a model using GAMS Direct"""

        options = self.options(options)
        config = self.config(config_options)

        gams = self._gams_module

        # What do I need to construct the GMO model?
        # list of all variables
        # list of all constraints, with linear coefficients, and nonlinear parts

        pass


GAMSDirect.solve.__doc__ = add_docstring_list(GAMSDirect.solve.__doc__, GAMSDirect.CONFIG)


if __name__ == "__main__":
    SolverFactory('new_gams_direct').solve(None)
