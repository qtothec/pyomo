import logging

from pyomo.core import Transformation, TransformationFactory

logger = logging.getLogger('pyomo.gdp.mpcc')


@TransformationFactory.register('gdp.mpcc', doc="Convert the GDP into MPCC")
class GDPtoMPCC(Transformation):
    """Convert a GDP into MPCC (Mathematical Program with Complementarity Constraints)

    First relax the model using a standard transformation to MINLP,
    then enforce integrality of the binaries using complementarity conditions.

    """
    pass
