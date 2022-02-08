from pyomo.environ import ConcreteModel, assert_optimal_termination, value, units as pyunits, check_optimal_termination
from idaes.core import FlowsheetBlock, MomentumBalanceType
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import fixed_variables_set, unfixed_variables_set, unused_variables_set
from pyomo.util.check_units import assert_units_consistent
from watertap.property_models.ion_DSPMDE_prop_pack import DSPMDEParameterBlock, ActivityCoefficientModel, DensityCalculation
from watertap.unit_models.nanofiltration_DSPMDE_0D import NanofiltrationDSPMDE0D, MassTransferCoefficient
from watertap.core.util.initialization import check_dof, check_solve
from watertap.core.util.infeasible import *
from idaes.core.util import get_solver
from watertap.core.util.infeasible import *
from idaes.core.util.model_statistics import *
import idaes.logger as idaeslog

import logging
from pyomo.util.infeasible import *


if __name__ == '__main__':
    solver = get_solver()

    m = ConcreteModel()



    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = DSPMDEParameterBlock(default={
        "solute_list": [
            "Ca_2+",
            "SO4_2-",
            "Mg_2+",
            "Na_+",
            "Cl_-",
            # 'X',
            # 'Y',
        ],
        "diffusivity_data": {
            ("Liq", "Ca_2+"): 0.792e-9,
            ("Liq", "SO4_2-"): 1.06e-9,
            ("Liq", "Mg_2+"): 0.706e-9,
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
            # ("Liq", "X"): 2.03e-9,
            # ("Liq", "Y"): 2.03e-9,

        },
        "mw_data": {
            "H2O": 18e-3,
            "Ca_2+": 40e-3,
            "Mg_2+": 24e-3,
            "SO4_2-": 96e-3,
            "Na_+": 23e-3,
            "Cl_-": 35e-3,
            # "X": 20e-3,
            # "Y": 20e-3,

        },
        "stokes_radius_data": {
            "Ca_2+": 0.309e-9,
            "Mg_2+": 0.347e-9,
            "SO4_2-": 0.230e-9,
            "Cl_-": 0.121e-9,
            "Na_+": 0.184e-9,
            # "X": 0.184e-9,
            # "Y": 0.184e-9,
        },
        "charge": {
            "Ca_2+": 2,
            "Mg_2+": 2,
            "SO4_2-": -2,
            "Na_+": 1,
            "Cl_-": -1,
            # "X": -1,
            # "Y": 1

        },
        "activity_coefficient_model": ActivityCoefficientModel.davies,
        "density_calculation": DensityCalculation.constant
    })

    m.fs.unit = NanofiltrationDSPMDE0D(default={"property_package": m.fs.properties,
                                                "mass_transfer_coefficient": MassTransferCoefficient.spiral_wound})
    mass_flow_in = 1 * pyunits.kg / pyunits.s
    feed_mass_frac = {
                      'Ca_2+': 382e-6,
                      'Mg_2+': 1394e-6,
                      'SO4_2-': 2136e-6,
                      'Cl_-': 20101.6e-6,
                      'Na_+': 11122e-6,
                      # 'Cl_-': 0.016924782608695656,
                      # 'X': 11122e-6,
                      # 'Y': 11122e-6,

                      }

    # Fix mole flow rates of each ion and water
    for ion, x in feed_mass_frac.items():
        mol_comp_flow = x * pyunits.kg / pyunits.kg * mass_flow_in / m.fs.unit.feed_side.properties_in[0].mw_comp[ion]
        m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', ion].fix(mol_comp_flow)
    H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
    H2O_mol_comp_flow = H2O_mass_frac * pyunits.kg / pyunits.kg * mass_flow_in / \
                        m.fs.unit.feed_side.properties_in[0].mw_comp['H2O']
    m.fs.unit.inlet.flow_mol_phase_comp[0, 'Liq', 'H2O'].fix(H2O_mol_comp_flow)

    # Use assert electroneutrality method from property model to ensure the ion concentrations provided
    # obey electroneutrality condition
    m.fs.unit.feed_side.properties_in[0].assert_electroneutrality(defined_state=True,
                                                                  adjust_by_ion='Cl_-',
                                                                  get_property='mass_frac_phase_comp')
    # m.fs.unit.feed_side.properties_in[0].mass_frac_phase_comp.display()

    # Fix other inlet state variables
    m.fs.unit.inlet.temperature[0].fix(298.15)
    m.fs.unit.inlet.pressure[0].fix(4e5)

    # Fix the membrane variables that are usually fixed for the DSPM-DE model
    m.fs.unit.radius_pore.fix(0.5e-9)
    m.fs.unit.membrane_thickness_effective.fix()
    m.fs.unit.membrane_charge_density.fix(-27)
    m.fs.unit.dielectric_constant_pore.fix(41.3)

    # Fix final permeate pressure to be ~atmospheric
    m.fs.unit.mixed_permeate[0].pressure.fix(101325)

    # Fix additional variables for calculating mass transfer coefficient with spiral wound correlation
    m.fs.unit.spacer_porosity.fix(1)
    m.fs.unit.spacer_mixing_efficiency.fix()
    m.fs.unit.spacer_mixing_length.fix()
    # m.fs.unit.length.fix(1)
    m.fs.unit.channel_height.fix(5e-4)
    # m.fs.unit.velocity[0, 0].fix(0.1)
    # m.fs.unit.recovery_vol_phase[0, 'Liq'].fix(0.5)
    # m.fs.unit.rejection_intrinsic_phase_comp[0, 'Liq', 'Ca_2+'].setlb(.2)



    # Set electroneutrality tolerance to 0 (value used in equality constraints for electroneutrality in unit model)
    m.fs.unit.tol_electroneutrality = 0

    iscale.calculate_scaling_factors(m.fs.unit)
    print('---------------- After calculate_scaling_factors---------------------------------------')
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]



    # checking state block
    assert_units_consistent(m)

    # check dof = 0
    check_dof(m, fail_flag=False)

    m.fs.unit.initialize()

    # Membrane area and recovery haven't typically been included in the literature for DSPM-DE,
    # but since we include them in our model to simulate/optimize at process level, I chose to fix those here


    # for con in m.fs.unit.component_data_objects(Constraint, descend_into=False):
    #     con.deactivate()
    # #
    # m.fs.unit.eq_water_flux.activate()
    m.fs.unit.eq_solute_solvent_flux.deactivate()
    # m.fs.unit.eq_solute_flux_concentration_polarization.activate()
    # m.fs.unit.eq_permeate_isothermal.activate()
    # m.fs.unit.eq_permeate_isothermal_mixed.activate()
    # m.fs.unit.eq_pressure_permeate_io.activate()
    # m.fs.unit.eq_mass_transfer_feed.activate()
    # m.fs.unit.eq_permeate_production.activate()
    # #
    # m.fs.unit.eq_solute_flux_pore_domain.deactivate()
    # m.fs.unit.eq_recovery_mol_phase_comp.activate()
    # m.fs.unit.eq_pore_isothermal.activate()
    # m.fs.unit.feed_side.eq_feed_interface_isothermal.activate()
    # m.fs.unit.feed_side.eq_feed_isothermal.activate()
    #
    # m.fs.unit.eq_recovery_vol_phase.activate()
    # m.fs.unit.eq_equal_flow_vol_pore_permeate.deactivate()
    # m.fs.unit.eq_equal_flow_vol_permeate.deactivate()
    #
    # # BRING IN NEW constraints for mass transfer coefficient
    # m.fs.unit.eq_Kf_comp.activate()
    # m.fs.unit.eq_N_Sc_comp.activate()
    # m.fs.unit.eq_N_Pe_comp.activate()
    # m.fs.unit.eq_velocity.activate()
    # m.fs.unit.eq_area.activate()
    #
    # m.fs.unit.eq_interfacial_partitioning_feed.activate()
    # m.fs.unit.eq_interfacial_partitioning_permeate.activate()
    # m.fs.unit.eq_electroneutrality_interface.activate()
    # m.fs.unit.eq_electroneutrality_pore.activate()
    # m.fs.unit.eq_electroneutrality_permeate.activate()
    # m.fs.unit.eq_rejection_intrinsic_phase_comp.activate()
    #
    # m.fs.unit.eq_equal_flowrate_pore_entrance_io.activate()
    # m.fs.unit.eq_pressure_pore_exit_io.activate()

    # Constraints that mess up the solve
    # m.fs.unit.eq_density_mixed_permeate.activate()
    # m.fs.unit.eq_equal_flowrate_pore_io.activate()
    # m.fs.unit.eq_equal_flowrate_feed_interface_in.activate()

    b = m.fs.unit

    # m.fs.unit.eq_equal_flow_vol_permeate.deactivate()
    # m.fs.unit.eq_equal_flow_vol_pore_permeate.deactivate()
    # m.fs.unit.recovery_vol_phase[0, 'Liq'].unfix()
    # m.fs.unit.area.unfix()

    # print('---------------- BEFORE AUTOMATE RESCALE---------------------------------------')
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
    print('\nNUMBER OF badly scaled variables:', len(list(iscale.badly_scaled_var_generator(m))))
    m.fs.unit._automate_rescale_variables(rescale_factor=1)
    # print('---------------- AFTER AUTOMATE RESCALE---------------------------------------')
    [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
    print('\nNUMBER OF badly scaled variables:', len(list(iscale.badly_scaled_var_generator(m))))


    # solver.options['max_iter'] = 0
    results = solver.solve(m, tee=True)
    if check_optimal_termination(results):
        b.report()
        print('SUCCESSFUL FIRST SOLVE!!!!!!!!!!!')
        # Check that electroneutrality is satisfied for feed outlet and mixed permeate- constraints that
        # are deactivated because they lead to failed solve
        b.feed_side.properties_out[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
        b.mixed_permeate[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
    else:
        print('FIRST SOLVE FAILED..............')
        print('Badly scaled vars after FIRST failed solve:')
        [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
        print(f'Number of badly scaled vars = {len(list(iscale.badly_scaled_var_generator(m)))}')
        m.fs.unit._automate_rescale_variables(rescale_factor=1)

        results = solver.solve(m, tee=True)
        if check_optimal_termination(results):
            b.report()
            print('SUCCESSFUL SECOND SOLVE!!!!!!!!!!!')
            # Check that electroneutrality is satisfied for feed outlet and mixed permeate- constraints that
            # are deactivated because they lead to failed solve
            b.feed_side.properties_out[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
            b.mixed_permeate[0].assert_electroneutrality(defined_state=False, tee=True, solve=False)
        else:
            print('Badly scaled vars after SECOND failed solve:')
            [print(i[0], i[1]) for i in iscale.badly_scaled_var_generator(m)]
            print(f'Number of badly scaled vars = {len(list(iscale.badly_scaled_var_generator(m)))}')
            print('SECOND SOLVE FAILED...............................')

    print(f'Degrees of freedom ={degrees_of_freedom(m)} ')

    print_infeasible_constraints(m, print_expression=True, print_variables=True, output_file='infeasible_w_kf.log')

