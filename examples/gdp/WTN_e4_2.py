# In terms of total flows and concentrations

from pyomo.environ import *
from pyomo.gdp import *

def build_water_treatment_network_model():
    """Build the water treatment network model"""
    m = ConcreteModel(name = "Water Treatment Network")

    
    """Set declarations"""
    m.in_flows = RangeSet(1, 3, doc="Inlet total flows", ordered=True)
    #Water is represented as fourth component, but really represents total flow
    m.comps = Set(initialize=['A', 'B', 'C', 'W'])
    m.mixers = RangeSet(1, 4, doc="Mixers", ordered=True)
    m.mixer_ins = RangeSet(1, 6, doc="Mixer_Ins", ordered=True)
    m.splitters = RangeSet(1, 6, doc="Splitters", ordered=True)
    m.splitter_outs = RangeSet(1, 4, doc="Splitter_Outs", ordered=True)
    m.tru = RangeSet(1, 3, doc="Treatment process units", ordered=True)
    

    """Parameter and initial point declarations"""

    #Inlet flow information
    #Component flow of water is just the same as the total flowrate
    #ppm for A,B and t/h for W
    inlet = {1: {'A':390, 'B':16780, 'C':25, 'W':13.1},
             2: {'A':10, 'B':110, 'C':100, 'W':32.7},
             3: {'A':25, 'B':40, 'C':35, 'W':56.5}}

    limits = {'A':2, 'B':2, 'C':5} # Discharge limits [=] ppm
    
    m.flow_total = sum(inlet[i]['W'] for i in m.in_flows)

    # equipment_info = {num: name, [removal ratio A, removal ratio B]}
    equipment_info = {1:['X', 99.9, 0.0, 0.0],
                      2:['XX', 90.0, 90.0, 97.0],
                      3:['XXX', 0.0, 95.0, 20.0]}
    
    @m.Param(m.tru, m.comps, doc="Equipment Removal Ratio for Each Component")
    def beta(m, equip, comp):
        if comp == 'A':
            return equipment_info[equip][1]/100
        elif comp == 'B':
            return equipment_info[equip][2]/100
        elif comp == 'C':
            return equipment_info[equip][3]/100
        else:
            return 0


    """Variable Declarations"""

    m.S_k = Var(m.splitters, m.splitter_outs, m.comps, domain=NonNegativeReals, doc="Splitter Effluent")
    m.M_k = Var(m.mixers, m.mixer_ins, m.comps, domain=NonNegativeReals, doc="Mixer Inlet")
    m.IPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Inlet")
    m.OPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Outlet")
    m.split = Var(m.splitters, m.splitter_outs, domain=NonNegativeReals, bounds=(0,1), 
                                      doc="Split fractions for splitter k into stream i")


    """Constraint definitions"""

    @m.Constraint(m.mixers, m.comps, doc="Flow Balance for mixer k")
    def mixer_balance(m, mixer, comp):
        if comp != 'W':
            if mixer < len(m.mixers):
                return sum(m.M_k[mixer,inlet,comp]*m.M_k[mixer,inlet,'W'] for inlet in m.mixer_ins) == m.IPU[mixer,comp]*m.IPU[mixer,'W']
            else: #last mixer has a different mass balance
                return sum(m.M_k[mixer,inlet,comp]*m.M_k[mixer,inlet,'W'] for inlet in m.mixer_ins) == limits[comp]*m.flow_total
        else:
            if mixer < len(m.mixers):
                return sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins) == m.IPU[mixer,comp]
            else:
                return sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins) == m.flow_total
    
    @m.Constraint(m.splitters, m.mixers, m.comps, doc="Splitter effluents are mixer inlets")
    def split_mix(m, splitter, mixer, comp):
            return m.M_k[mixer,splitter,comp] == m.S_k[splitter,mixer,comp]

    @m.Constraint(m.splitters, doc="Split Fraction Balance for splitter k")
    def split_fraction_balance(m, splitter):
        return 1 == sum(m.split[splitter,outlet] for outlet in m.splitter_outs)

    @m.Constraint(m.splitters, m.splitter_outs, m.comps, doc="Component Split Balance for splitter k")
    def splitter_balance(m, splitter, outlet, comp):
        if comp != 'W':
            if splitter <= len(m.in_flows):
                return inlet[splitter][comp] == m.S_k[splitter,outlet,comp]
            else:
                return m.OPU[splitter-len(m.in_flows),comp] == m.S_k[splitter,outlet,comp]
        else:
            if splitter <= len(m.in_flows): #based on number of inlet streams
                return inlet[splitter][comp] == sum(m.S_k[splitter,outlet,comp] for outlet in m.splitter_outs)
            else:
                return m.OPU[splitter-len(m.in_flows),comp] == sum(m.S_k[splitter,outlet,comp] for outlet in m.splitter_outs)

    @m.Constraint(m.tru, m.comps, doc="Component Removal for Treatment Unit k")
    def component_removal(m,equip,comp):
        return m.OPU[equip,comp] == (1-m.beta[equip,comp])*m.IPU[equip,comp]
        
    
    """Objective function definition"""
    
    m.minCost = Objective(expr=sum(m.IPU[t,'W'] for t in m.tru), doc="Minimize waste stream processing cost")

    return m


model = build_water_treatment_network_model()

opt = SolverFactory('gams')
results = opt.solve (model, tee=True, solver='baron')

print(results)
model.pprint()
