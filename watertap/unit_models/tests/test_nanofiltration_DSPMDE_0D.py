###############################################################################
# WaterTAP Copyright (c) 2021, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National
# Laboratory, National Renewable Energy Laboratory, and National Energy
# Technology Laboratory (subject to receipt of any required approvals from
# the U.S. Dept. of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#
###############################################################################

import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           Var,
                           units as pyunits,
                           assert_optimal_termination)
from pyomo.network import Port
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
import watertap.property_models.ion_DSPMDE_prop_pack as props
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D
from watertap.core.util.initialization import check_dof

from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (calculate_scaling_factors,
                                     unscaled_variables_generator,
                                     unscaled_constraints_generator,
                                     badly_scaled_var_generator)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

class TestNanoFiltration():
    @pytest.fixture(scope="class")
    def NF_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = props.DSPMDEParameterBlock(default={
            "solute_list": ["Ca_2+", "SO4_2-", "Na_+", "Cl_-", "Mg_2+"],
            "diffusivity_data": {("Liq", "Ca_2+"): 0.792e-9,
                                 ("Liq", "SO4_2-"): 1.06e-9,
                                 ("Liq", "Na_+"): 1.33e-9,
                                 ("Liq", "Cl_-"): 2.03e-9,
                                 ("Liq", "Mg_2+"): 0.706e-9},
            "mw_data": {"H2O": 18e-3,
                        "Na_+": 23e-3,
                        "Ca_2+": 40e-3,
                        "Mg_2+": 24e-3,
                        "Cl_-": 35e-3,
                        "SO4_2-": 96e-3},
            "stokes_radius_data": {"Na_+": 0.184e-9,
                                   "Ca_2+": 0.309e-9,
                                   "Mg_2+": 0.347e-9,
                                   "Cl_-": 0.121e-9,
                                   "SO4_2-": 0.230e-9},
            "density_data": {"H2O": 1000,
                             "Na_+": 968,
                             "Ca_2+": 1550,
                             "Mg_2+": 1738,
                             "Cl_-": 3214,
                             "SO4_2-": 2553},
            "charge": {"Na_+": 1,
                       "Ca_2+": 2,
                       "Mg_2+": 2,
                       "Cl_-": -1,
                       "SO4_2-": -2},
            })

        m.fs.unit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties})

        return m

    @pytest.mark.component
    def test_property_ions(self, NF_frame):
        m = NF_frame
        m.fs.stream = m.fs.properties.build_state_block([0], default={'defined_state': True})


        mass_flow_in = 1 * pyunits.kg/pyunits.s
        feed_mass_frac = {'Na_+': 11122e-6,
                          'Ca_2+': 382e-6,
                          'Mg_2+': 1394e-6,
                          'SO4_2-': 2136e-6,
                          'Cl_-': 20300e-6}
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = x * pyunits.kg/pyunits.kg * mass_flow_in / m.fs.stream[0].mw_comp[ion]
            m.fs.stream[0].flow_mol_phase_comp['Liq', ion].fix(mol_comp_flow)

        m.fs.stream[0].temperature.fix(298.15)
        m.fs.stream[0].pressure.fix(101325)
        m.fs.stream[0].assert_electroneutrality(tol=1e-2)

    @pytest.mark.unit
    def test_config(self, NF_frame):
        m = NF_frame
        mass_flow_in = 1 * pyunits.kg/pyunits.s
        feed_mass_frac = {'Na_+': 11122e-6,
                          'Ca_2+': 382e-6,
                          'Mg_2+': 1394e-6,
                          'SO4_2-': 2136e-6,
                          'Cl_-': 20300e-6}
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = x * pyunits.kg/pyunits.kg * mass_flow_in / m.fs.stream[0].mw_comp[ion]
            m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', ion].fix(mol_comp_flow)
            m.fs.unit.permeate.flow_mol_phase_comp[0, 'Liq', ion].fix(0)#TODO: remove later- temporary to eliminate dof


        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in / \
                            m.fs.unit.feed_side.properties_in[0].mw_comp['H2O']

        m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(H2O_mol_comp_flow)
        m.fs.unit.inlet.temperature[0].fix(298.15)
        m.fs.unit.inlet.pressure[0].fix(101325)
        m.fs.unit.radius_pore.fix(0.5e-9)
        m.fs.unit.permeate.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(0) #TODO: remove later- temporary to eliminate dof

        # check unit config arguments
        assert len(m.fs.unit.config) == 7

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert m.fs.unit.config.material_balance_type == \
               MaterialBalanceType.useDefault
        assert m.fs.unit.config.momentum_balance_type == \
               MomentumBalanceType.pressureTotal
        assert not m.fs.unit.config.has_pressure_change
        assert m.fs.unit.config.property_package is \
               m.fs.properties

    # @pytest.mark.unit
    # def test_build(self, NF_frame):
    #     m = NF_frame
    #
    #     # test ports and variables
    #     port_lst = ['inlet', 'retentate', 'permeate']
    #     port_vars_lst = ['flow_mass_phase_comp', 'pressure', 'temperature']
    #     for port_str in port_lst:
    #         assert hasattr(m.fs.unit, port_str)
    #         port = getattr(m.fs.unit, port_str)
    #         assert len(port.vars) == 3
    #         assert isinstance(port, Port)
    #         for var_str in port_vars_lst:
    #             assert hasattr(port, var_str)
    #             var = getattr(port, var_str)
    #             assert isinstance(var, Var)
    #
    #     # test unit objects (including parameters, variables, and constraints)
    #     unit_objs_lst = ['A_comp', 'B_comp', 'sigma', 'dens_solvent',
    #                      'flux_mass_phase_comp_in', 'flux_mass_phase_comp_out',
    #                      'avg_conc_mass_phase_comp_in', 'avg_conc_mass_phase_comp_out', 'area',
    #                      'deltaP', 'mass_transfer_phase_comp', 'flux_mass_phase_comp_avg',
    #                      'eq_mass_transfer_term', 'eq_permeate_production',
    #                      'eq_flux_in', 'eq_flux_out',
    #                      'eq_avg_conc_in', 'eq_avg_conc_out',
    #                      'eq_connect_mass_transfer', 'eq_connect_enthalpy_transfer',
    #                      'eq_permeate_isothermal']
    #     for obj_str in unit_objs_lst:
    #         assert hasattr(m.fs.unit, obj_str)
    #
    #
    #     # test state block objects
    #     cv_name = 'feed_side'
    #     cv_stateblock_lst = ['properties_in', 'properties_out']
    #     stateblock_objs_lst = \
    #         ['flow_mass_phase_comp', 'pressure', 'temperature', 'pressure_osm',
    #          'osm_coeff', 'mass_frac_phase_comp', 'conc_mass_phase_comp',
    #          'dens_mass_phase', 'enth_mass_phase',
    #          'eq_pressure_osm', 'eq_osm_coeff', 'eq_mass_frac_phase_comp',
    #          'eq_conc_mass_phase_comp', 'eq_dens_mass_phase', 'eq_enth_mass_phase'
    #          ]
    #     # control volume
    #     assert hasattr(m.fs.unit, cv_name)
    #     cv_blk = getattr(m.fs.unit, cv_name)
    #     for blk_str in cv_stateblock_lst:
    #         assert hasattr(cv_blk, blk_str)
    #         blk = getattr(cv_blk, blk_str)
    #         for obj_str in stateblock_objs_lst:
    #             assert hasattr(blk[0], obj_str)
    #     # permeate
    #     assert hasattr(m.fs.unit, 'properties_permeate')
    #     blk = getattr(m.fs.unit, 'properties_permeate')
    #     for var_str in stateblock_objs_lst:
    #         assert hasattr(blk[0], var_str)
    #
    #     # test statistics
    #     assert number_variables(m) == 73
    #     assert number_total_constraints(m) == 45
    #     assert number_unused_variables(m) == 7  # vars from property package parameters
    #
    @pytest.mark.unit
    def test_dof(self, NF_frame):
        m = NF_frame
        check_dof(m, fail_flag=True)

    @pytest.mark.unit
    def test_calculate_scaling(self, NF_frame):
        m = NF_frame
        # calculate_scaling_factors(m.fs.unit)

        # check that all variables have scaling factors
        # unscaled_var_list = list(unscaled_variables_generator(m))
        # assert len(unscaled_var_list) == 0
        # # check that all constraints have been scaled
        # unscaled_constraint_list = list(unscaled_constraints_generator(m))
        # assert len(unscaled_constraint_list) == 0

    @pytest.mark.component
    def test_initialize(self, NF_frame):
        m= NF_frame
        # initialization_tester(m)
        m.fs.unit.initialize()

    # @pytest.mark.component
    # def test_var_scaling(self, NF_frame):
    #     m = NF_frame
    #     badly_scaled_var_lst = list(badly_scaled_var_generator(m))
    #     [print(i,j) for i, j in badly_scaled_var_lst]
    #     assert badly_scaled_var_lst == []
    #
    @pytest.mark.component
    def test_solve(self, NF_frame):
        m = NF_frame
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)
    #
    # @pytest.mark.component
    # def test_conservation(self, NF_frame):
    #     m = NF_frame
    #     b = m.fs.unit
    #     comp_lst = ['NaCl', 'H2O']
    #
    #     flow_mass_inlet = sum(
    #         b.feed_side.properties_in[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #     flow_mass_retentate = sum(
    #         b.feed_side.properties_out[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #     flow_mass_permeate = sum(
    #         b.properties_permeate[0].flow_mass_phase_comp['Liq', j] for j in comp_lst)
    #
    #     assert (abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate
    #                       )) <= 1e-6)
    #
    #     assert (abs(value(
    #         flow_mass_inlet * b.feed_side.properties_in[0].enth_mass_phase['Liq']
    #         - flow_mass_retentate * b.feed_side.properties_out[0].enth_mass_phase['Liq']
    #         - flow_mass_permeate * b.properties_permeate[0].enth_mass_phase['Liq']
    #     )) <= 1e-6)
    #
    # @pytest.mark.component
    # def test_solution(self, NF_frame):
    #     m = NF_frame
    #     assert (pytest.approx(1.079e-2, rel=1e-3) ==
    #             value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'H2O']))
    #     assert (pytest.approx(3.435e-4, rel=1e-3) ==
    #             value(m.fs.unit.flux_mass_phase_comp_avg[0, 'Liq', 'NaCl']))
    #     assert (pytest.approx(0.5396, rel=1e-3) ==
    #             value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'H2O']))
    #     assert (pytest.approx(1.717e-2, rel=1e-3) ==
    #             value(m.fs.unit.properties_permeate[0].flow_mass_phase_comp['Liq', 'NaCl']))
