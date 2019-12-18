from pyutilib.misc.config import ConfigValue

from pyomo.common.config import add_docstring_list
from pyomo.solver.base import OptSolver


class GAMSDirect(OptSolver):
    """Direct interface to GAMS"""
    CONFIG = OptSolver.CONFIG()
    MAPPED_OPTIONS = OptSolver.MAPPED_OPTIONS()

    CONFIG.declare('symbolic_solver_labels', ConfigValue(
        default=False, domain=bool,
        doc="If True, the GAMS variable and constraint names will "
        "closely match those of the pyomo variables and constraints."
    ))
    CONFIG.declare('stream_solver', ConfigValue(
        default=False, domain=bool,
        doc="If True, show the GAMS output"
    ))
    CONFIG.declare('load_solutions', ConfigValue(
        default=True, domain=bool,
        doc="If True, load the solution back into the Pyomo model"
    ))

    __doc__ = add_docstring_list(__doc__, CONFIG)

    def __init__(self):
        super(GAMSDirect, self).__init__()

        try:
            import gams
        finally:
            pass
