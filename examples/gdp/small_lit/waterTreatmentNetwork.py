'''
Problem from Example 3.5.1: Water Treatment Network
(Generalized Disjunctive Programming: A Framework for Formulation and
 Alternative Algorithms for MINLP Optimization - Grossmann and Ruiz, 2012)

This example corresponds to a synthesis problem of a distributed
wastewater multicomponent network. Given a set of process liquid
streams with known composition, a set of technologies for the removal
of pollutants, and a set of mixers and splitters, the objective is to
find the interconnections of the technologies and their flowrates to 
meet the specified discharge composition of pollutant at minimum
total cost. Discrete choices involve deciding what equipment to use for 
each treatment unit.

Link to Article: 
https://pdfs.semanticscholar.org/d533/63bc770f93f876277136a5c6a9fba12a27e6.pdf

Original Problem from Galan and Grossmann, 1998, Example 1

Full GDP model and explanation of model for wastewater treatment network from
Lee and Grossmann (2003)
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
    #Water is represented as fourth component, but really represents total flow
    m.comps = Set(initialize=['A', 'B', 'C', 'W'])
    m.mixers = RangeSet(1, 4, doc="Mixers", ordered=True)
    m.mixer_ins = RangeSet(1,6, doc="Mixer_Ins", ordered=True)
    m.splitters = RangeSet(1, 6, doc="Splitters", ordered=True)
    m.splitter_outs = RangeSet(1, 4, doc="Splitter_Outs", ordered=True)
    m.tru = RangeSet(1, 3, doc="Treatment process units", ordered=True)
    m.allEquipment = RangeSet(1, 9, doc="All Equipment", ordered=True)
    

    """Parameter and initial point declarations"""
    
    #Inlet flow information
    in_flows = {1:20, 2:15, 3:5} # t/h
    #Component flow of water is just the same as the total flowrate
    in_concs = {1: {'A':1100, 'B':300, 'C':400, 'W':1}, #ppm
                2: {'A':300, 'B':700, 'C':1500, 'W':1}, 
                3: {'A':500, 'B':1000, 'C':600, 'W':1}}

    @m.Param(m.in_flows, m.comps, doc="Inlet Component Flows [=] t*ppm/h")
    def in_comp_flow(m, flow, comp):
        return in_flows[flow] * (in_concs[flow][comp])
    
    limits = {'A':100,'B':100,'C':100, 'W':1} # Discharge limits [=] ppm
    
    m.out_flow_total = sum(in_flows[i] for i in m.in_flows)
    @m.Param(m.comps, doc="Outlet Component Flows [=] t*ppm/h")
    def out_comp_flow(m, comp):
        return m.out_flow_total * limits[comp]

    # equipment_info = {num: name, unit, removal ratio A, removal ratio B, removal ratio C, investment alpha, operating gamma}
    equipment_info = {1:['EA', 1, 90.0,  0.0, 40.0, 3480,      0],
                      2:['EB', 1, 50.0, 70.0,  0.0,  469,     10],
                      3:['EC', 1,  0.0, 80.0,  0.0,   26,      1],
                      4:['ED', 2,  0.0, 90.0,  0.0,  726, 0.0089],
                      5:['EE', 2,  0.0, 99.0,  0.0, 1260,  0.018],
                      6:['EF', 2, 50.0, 99.0, 80.0, 5000,    5.8],
                      7:['EG', 3, 80.0,  0.0, 60.0,  320,      6],
                      8:['EH', 3,  0.0,  0.0, 80.0,   58,     15],
                      9:['EI', 3,  0.0,  0.0, 40.0,   10,      1]}

    @m.Param(m.allEquipment, m.comps, doc="Equipment Removal Ratio for Each Component")
    def beta(m, equip, comp):
        if comp == 'A':
            return equipment_info[equip][2]/100
        elif comp == 'B':
            return equipment_info[equip][3]/100
        elif comp == 'C':
            return equipment_info[equip][4]/100
        else:
            return 0

    @m.Param(m.allEquipment, doc="Equipment Investment Cost")
    def alpha(m, equip):
        return equipment_info[equip][5]

    @m.Param(m.allEquipment, doc="Equipment Operating Cost")
    def gamma(m, equip):
        return equipment_info[equip][6]


    """Variable Declarations"""

    m.S_k = Var(m.splitters, m.splitter_outs, m.comps, domain=NonNegativeReals, doc="Splitter Effluent Streams")
    m.M_k = Var(m.mixers, m.mixer_ins, m.comps, domain=NonNegativeReals, doc="Mixer Inlet Streams")
    m.IPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Inlet Streams")
    m.OPU = Var(m.tru, m.comps, domain=NonNegativeReals, doc="TRU Outlet Streams")
    m.split = Var(m.splitters, m.splitter_outs, domain=NonNegativeReals, bounds=(0,1), 
                                      doc="Split fractions for splitter k into stream i")
    m.CP_k = Var(m.tru, domain=NonNegativeReals, doc="Cost of equipment h chosen for treatment unit k")
    
    m.S_k[1,1,'W'].fix(0)
    m.S_k[1,4,'W'].fix(0)
    m.S_k[2,1,'W'].fix(0)
    m.S_k[2,2,'W'].fix(15)
    m.S_k[2,3,'W'].fix(0)
    m.S_k[2,4,'W'].fix(0)
    m.S_k[3,1,'W'].fix(0)
    m.S_k[3,2,'W'].fix(5)
    m.S_k[3,3,'W'].fix(0)
    m.S_k[3,4,'W'].fix(0)
    m.S_k[4,1,'W'].fix(0)
    m.S_k[4,2,'W'].fix(0)
    m.S_k[5,2,'W'].fix(0)
    m.S_k[5,3,'W'].fix(0)
    m.S_k[5,4,'W'].fix(0)
    m.S_k[6,1,'W'].fix(0)
    m.S_k[6,2,'W'].fix(0)
    m.S_k[6,3,'W'].fix(0)
    
    """Constraint definitions"""

    @m.Constraint(m.mixers, m.comps, doc="Flow Balance for mixer k")
    def mixer_balance(m, mixer, comp):
        if mixer < len(m.mixers):
            return m.IPU[mixer,comp] == sum(m.M_k[mixer,inlet,comp] for inlet in m.mixer_ins)
        else: #last mixer has a different mass balance
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
                           
    @m.Constraint(m.tru, doc="Total Flow Balance for Treatment Unit k")
    def total_tru_flow_balance(m, equip):
        return m.IPU[equip,'W'] == m.OPU[equip,'W']
        
    """Disjunctions"""

    #Different types of equipment can be used for each treatment unit
    #Technologies have constant removal ratios for each pollutant

    @m.Disjunct(m.allEquipment)
    def equipment_disjuncts(disj,equip):

        @disj.Constraint(m.comps)
        def component_removal(disj,comp):
            return m.OPU[((equip-1)//3)+1,comp] == (1-m.beta[equip,comp])*m.IPU[((equip-1)//3)+1,comp]
        
        F = m.OPU[((equip-1)//3)+1,'W']

        disj.cost = Constraint(expr=m.CP_k[((equip-1)//3)+1] == (m.alpha[equip]*(F**0.7) + m.gamma[equip]*F))

    m.TX = Disjunction(expr=[m.equipment_disjuncts[1], m.equipment_disjuncts[2], m.equipment_disjuncts[3]], 
                           doc="Treatment Unit 1")

    m.TXX = Disjunction(expr=[m.equipment_disjuncts[4], m.equipment_disjuncts[5], m.equipment_disjuncts[6]], 
                           doc="Treatment Unit 2")

    m.TXXX = Disjunction(expr=[m.equipment_disjuncts[7], m.equipment_disjuncts[8], m.equipment_disjuncts[9]], 
                           doc="Treatment Unit 3")


    """Objective function definition"""
    
    m.minCost = Objective(expr=sum(m.CP_k[i] for i in m.tru), doc="Minimize waste stream processing cost")

    return m


model = build_water_treatment_network_model()

TransformationFactory('gdp.bigm').apply_to(model,bigM=1e6)

opt = SolverFactory('gams')

results = opt.solve (model, tee=True, solver='baron', add_options=['option reslim=120;'])

print(results)

model.pprint()


