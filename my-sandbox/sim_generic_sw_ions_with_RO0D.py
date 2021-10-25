###############################################################################
# ProteusLib Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/nawi-hub/proteuslib/"
#
###############################################################################
import pyomo.opt
from pyomo.environ import ConcreteModel, value, Constraint, units as pyunits, Var, Objective, SolverFactory
from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from pyomo.util.check_units import assert_units_consistent
from watertap.util.initialization import check_dof
from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from watertap.property_models.seawater_ion_generic import configuration
from watertap.unit_models.reverse_osmosis_0D import ReverseOsmosis0D, ConcentrationPolarizationType, MassTransferCoefficient, PressureChangeType, EnergyBalanceType
from idaes.core.util import get_solver
from watertap.util.initialization import *
from idaes.core.util.model_diagnostics import DegeneracyHunter

solver = get_solver()
# -----------------------------------------------------------------------------

m = ConcreteModel()

m.fs = FlowsheetBlock(default={"dynamic": False})
m.fs.properties = GenericParameterBlock(default=configuration)
m.fs.unit = ReverseOsmosis0D(default={
        "property_package": m.fs.properties,
        "energy_balance_type": EnergyBalanceType.none,
        "has_pressure_change": False,
        "concentration_polarization_type": ConcentrationPolarizationType.none,
        "mass_transfer_coefficient": MassTransferCoefficient.none,
        "pressure_change_type": PressureChangeType.fixed_per_stage})

# specify
#todo: checking electroneutrality at inlet --> sum(charge*conc_mol_phase_comp) = -0.1082
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'Na_+'].fix(0.008845)
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'Ca_2+'].fix(0.000174)
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'Mg_2+'].fix(0.001049)
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'SO4_2-'].fix(0.000407)
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'Cl_-'].fix(0.010479)
m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0.979046)
m.fs.unit.inlet.temperature.fix(273.15 + 25)
m.fs.unit.inlet.pressure.fix(50e5)

m.fs.unit.recovery_vol_phase.fix(0.4)
m.fs.unit.permeate_side.properties_mixed[0].pressure.fix(101325)
# m.fs.unit.permeate_side.properties_mixed[0].temperature.fix(273.15 + 25)

B= 3.5e-8
m.fs.unit.A_comp[0,'H2O'].fix(4.2e-12)
m.fs.unit.B_comp[0, 'Na_+'].set_value(B)
m.fs.unit.B_comp[0, 'Ca_2+'].set_value(B)
m.fs.unit.B_comp[0, 'Mg_2+'].set_value(B)
m.fs.unit.B_comp[0, 'SO4_2-'].set_value(B)
m.fs.unit.B_comp[0, 'Cl_-'].set_value(B)  # guess, but electroneutrality enforced below
m.fs.unit.rejection_phase_comp[0, 'Liq','Na_+'].fix(.98)
m.fs.unit.rejection_phase_comp[0, 'Liq','Cl_-'].set_value(.98)
m.fs.unit.rejection_phase_comp[0, 'Liq','Ca_2+'].fix(.999)
m.fs.unit.rejection_phase_comp[0, 'Liq','Mg_2+'].fix(.999)
m.fs.unit.rejection_phase_comp[0, 'Liq', 'SO4_2-'].fix(.999)
# m.fs.unit.rejection_phase_comp[0, 'Liq','Cl_-'].setlb(.97)
# m.fs.unit.deltaP.fix(0)

m.fs.unit.area.fix(50)
# m.fs.unit.area.setub(100)
charge_comp = {'Na_+': 1, 'Ca_2+': 2, 'Mg_2+': 2, 'SO4_2-': -2, 'Cl_-': -1}
m.fs.unit.eq_electroneutrality = Constraint(
    expr=0 ==
         sum(charge_comp[j]
             * m.fs.unit.feed_side.properties_out[0].conc_mol_phase_comp['Liq', j]
             for j in charge_comp))
# iscale.constraint_scaling_transform(m.fs.unit.eq_electroneutrality, 1)
# m.fs.stream[0].flow_mol_phase_comp
# m.fs.stream[0].mw_comp
# m.fs.stream[0].flow_mass_phase_comp
# m.fs.stream[0].pressure_osm_phase
# m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'Na_+'].fix(0.008845)
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'Ca_2+'].fix(0.000174)
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'Mg_2+'].fix(0.001049)
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'SO4_2-'].fix(0.000407)
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'Cl_-'].fix(0.010479)
# m.fs.stream[0].flow_mol_phase_comp['Liq', 'H2O'].fix(0.979046)
# m.fs.stream[0].temperature.fix(273.15 + 25)
# m.fs.stream[0].pressure.fix(50e5)
# # m.fs.stream[0].molality_phase_comp_apparent
# m.fs.stream[0].mass_frac_phase_comp

# m.fs.unit.feed_side.properties_in[0].mw_comp['H2O']

# # scaling
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1, index=('Liq', 'H2O'))
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Na'))
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e4, index=('Liq', 'Ca'))
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'Mg'))
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e3, index=('Liq', 'SO4'))
# m.fs.properties.set_default_scaling('flow_mass_phase_comp', 1e2, index=('Liq', 'Cl'))

m.fs.properties.set_default_scaling('pressure_osm_phase', 1e-5, index=('Liq'))
# Set scaling below based on printed results from Ben's perturbation_helper (to eliminate IPOPT induced changes)
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'NaCl'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'CaCl2'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'CaSO4'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'Na2SO4'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'MgCl2'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp_apparent', 1e-5, index=('Liq', 'MgSO4'))
m.fs.properties.set_default_scaling('mole_frac_phase_comp', 1e2, index=('Liq', 'H2O'))


# m.fs.properties.set_default_scaling('mole_frac_phase_comp', 1e4, index=('Liq', 'Cl'))
# m.fs.unit.osm_press_min=Constraint(expr=m.fs.unit.permeate_side.properties_mixed[0].pressure_osm_phase['Liq']
#                                         >=1e-2 * pyunits.kg / pyunits.m / pyunits.s**2)

iscale.calculate_scaling_factors(m)
#
# # checking state block
# assert_units_consistent(m.fs)

check_dof(m, fail_flag=True)


#
# # # initialize
# m.fs.properties.initialize()

m.fs.unit.initialize(fail_on_warning=False, initialize_guess={'cp_modulus': 1})
# print_initialization_perturbation(m.fs, 1e-8, 1e-2, 1e-8, True)

# # # # solve

solver.options = {'bound_push': 1e-8, 'halt_on_ampl_error': 'yes'}
# Specifying an iteration limit of 0 allows us to inspect the initial point
# solver.options['max_iter'] = 0

# "Solving" the model with an iteration limit of 0 load the initial point and applies
# any preprocessors (e.g., enforces bounds)
# m.fs.objective = Objective(expr=-m.fs.unit.recovery_vol_phase[0, "Liq"])

solver.solve(m, symbolic_solver_labels=True, tee=True)

# dh = DegeneracyHunter(m, solver=pyomo.opt.SolverFactory('ipopt'))

# solver.solve(m, tee=True, symbolic_solver_labels=True)

# m.fs.unit.initialize(fail_on_warning=False)

