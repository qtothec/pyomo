'''
Problem from Example 1 (Optimal Design of Distributed Wastewater Treatment Networks - Galan and Grossmann, 1998)

Link to Article: https://pubs.acs.org/doi/10.1021/ie980133h 

This example corresponds to a synthesis problem of a distributed
wastewater multicomponent network. Given a set of process liquid
streams with known composition, a set of technologies for the removal
of pollutants, and a set of mixers and splitters, the objective is to
find the interconnections of the technologies and their flowrates to 
meet the specified discharge composition of pollutant at minimum
total cost.

In this simplified example of the wastewater treatment network problem, costs are optimized by minimizing the flow rate treated and treatment unit equipment is pre-selected (no discrete choices involved).

Full GDP model and explanation of model for the wastewater treatment network problem can be found under example 4 of Lee and Grossmann (2003): 

https://ac-els-cdn-com.proxy.library.cmu.edu/S009813540300098X/1-s2.0-
S009813540300098X-main.pdf?_tid=b8fa9ebf-0035-4764-8cbd-aff04fe33274&acdnat=
1544402786_0af0a426c9d1438838d007e0cd023399
'''

from pyomo.environ import *
from pyomo.gdp import *

def build_water_treatment_network_model():
    """Build the water treatment network model"""
    m = ConcreteModel(name = "Water Treatment Network")

    
    """Set declarations"""
    m.in_flows = RangeSet(1, 3, doc="Inlet total flows", ordered=True)
    #Water is represented as third component, but really represents total flow
    m.comps = Set(initialize=['A', 'B', 'C', 'W'])
    m.mixers = RangeSet(1, 4, doc="Mixers", ordered=True)
    m.mixer_ins = RangeSet(1, 6, doc="Mixer_Ins", ordered=True)
    m.splitters = RangeSet(1, 6, doc="Splitters", ordered=True)
    m.splitter_outs = RangeSet(1, 4, doc="Splitter_Outs", ordered=True)
    m.tru = RangeSet(1, 3, doc="Treatment process units", ordered=True)
    

    """Parameter and initial point declarations"""

    #Inlet flow information
    in_flows = {1:13.1, 2:32.7, 3:56.5} # t/h
    #Component flow of water is just the same as the total flowrate
    in_concs = {1: {'A':390, 'B':16780, 'C':25, 'W':1}, #ppm
                2: {'A':10, 'B':110, 'C':100, 'W':1},
                3: {'A':25, 'B':40, 'C':35, 'W':1}}

    @m.Param(m.in_flows, m.comps, doc="Inlet Component Flows [=] t*ppm/h")
    def in_comp_flow(m, flow, comp):
        return in_flows[flow] * (in_concs[flow][comp])

    limits = {'A':2, 'B':2, 'C':5, 'W':1} # Discharge limits [=] ppm
    
    m.out_flow_total = sum(in_flows[i] for i in m.in_flows)
    @m.Param(m.comps, doc="Outlet Component Flows [=] t*ppm/h")
    def out_comp_flow(m, comp):
        return m.out_flow_total * limits[comp]

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

    m.S_k = Var(m.splitters, m.splitter_outs, m.comps, domain=NonNegativeReals, doc="Splitter Effluent Streams")
    m.M_k = Var(m.mixers, m.mixer_ins, m.comps, domain=NonNegativeReals, doc="Mixer Inlet Streams")
    m.IPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Inlet Streams")
    m.OPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Outlet Streams")
    m.split = Var(m.splitters, m.splitter_outs, domain=NonNegativeReals, bounds=(0,1), 
                                      doc="Split fractions for splitter k into stream i")


    """Constraint definitions"""

    @m.Constraint(m.mixers, m.comps, doc="Flow Balance for mixer k")
    def mixer_balance(m, mixer, comp):
        if mixer < len(m.mixers):
            return m.IPU[mixer,comp] == sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins)
        else: 
            #last mixer has a different mass balance
            if comp == 'W':
                return m.out_comp_flow[comp] == sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins)
            else:  
                return m.out_comp_flow[comp] >= sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins)
    
    @m.Constraint(m.splitters, m.mixers, m.comps, doc="Splitter effluents are mixer inlets")
    def split_mix(m, splitter, mixer, comp):
        return m.M_k[mixer,splitter,comp] == m.S_k[splitter,mixer,comp]
    
    @m.Constraint(m.splitters, m.comps, doc="Flow Balance for splitter k")
    def splitter_balance(m, splitter, comp):
        if splitter <= len(m.in_flows): #based on number of inlet streams
            return m.in_comp_flow[splitter,comp] == sum(m.S_k[splitter,outlet,comp] for outlet in m.splitter_outs)
        else:
            return m.OPU[splitter-len(m.in_flows),comp] == sum(m.S_k[splitter,outlet,comp] for outlet in m.splitter_outs)

    @m.Constraint(m.splitters, doc="Split Fraction Balance for splitter k")
    def split_fraction_balance(m, splitter):
        return 1 == sum(m.split[splitter,outlet] for outlet in m.splitter_outs)

    @m.Constraint(m.splitters, m.splitter_outs, m.comps, doc="Flow Split Balance for splitter k")
    def flow_split_balance(m, splitter, outlet, comp):
        if splitter <= len(m.in_flows):
            return m.S_k[splitter,outlet,comp] == m.split[splitter,outlet] * m.in_comp_flow[splitter,comp]
        else:
            return m.S_k[splitter,outlet,comp] == m.split[splitter,outlet] * m.OPU[splitter-len(m.in_flows),comp]

    @m.Constraint(m.tru, m.comps, doc="Component Removal for Treatment Unit k")
    def component_removal(m,equip,comp):
        return m.OPU[equip,comp] == (1-m.beta[equip,comp])*m.IPU[equip,comp]

    @m.Constraint(m.tru, doc="Total Flow Balance for Treatment Unit k")
    def total_tru_flow_balance(m, equip):
        return m.IPU[equip,'W'] == m.OPU[equip,'W']
    
    """Objective function definition"""
    
    m.minCost = Objective(expr=sum(m.IPU[t,'W'] for t in m.tru), doc="Minimize waste stream processing cost")

    return m


model = build_water_treatment_network_model()

opt = SolverFactory('gams')
results = opt.solve(model, tee=True, solver='baron')

print(results)
model.pprint()
