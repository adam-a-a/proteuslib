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

# Importing testing libraries
import pytest

# Importing the object for units from pyomo
from pyomo.environ import units as pyunits, Reals, Param

# Imports from idaes core
from idaes.core import AqueousPhase
from idaes.core.components import Solvent, Solute, Cation, Anion, Apparent
from idaes.core.phases import PhaseType as PT

# Imports from idaes generic models
import idaes.generic_models.properties.core.pure.Perrys as Perrys
from idaes.generic_models.properties.core.pure.electrolyte import relative_permittivity_constant
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.eos.enrtl_reference_states import Unsymmetric

# Importing the enum for concentration unit basis used in the 'get_concentration_term' function
from idaes.generic_models.properties.core.generic.generic_reaction import ConcentrationForm

# Import the object/function for heat of reaction
from idaes.generic_models.properties.core.reactions.dh_rxn import constant_dh_rxn

# Import safe log power law equation
from idaes.generic_models.properties.core.reactions.equilibrium_forms import log_power_law_equil

# Import built-in van't Hoff function
from idaes.generic_models.properties.core.reactions.equilibrium_constant import van_t_hoff

# Import specific pyomo objects
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Suffix)

# Import the scaling methods
from idaes.core.util import scaling as iscale

# Import pyomo methods to check the system units
from pyomo.util.check_units import assert_units_consistent

# Import idaes methods to check the model during construction
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)

# Import the idaes objects for Generic Properties and Reactions
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock, StateIndex)
from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)

# Import the idaes object for the EquilibriumReactor unit model
from idaes.generic_models.unit_models.equilibrium_reactor import EquilibriumReactor

# Import the core idaes objects for Flowsheets and types of balances
from idaes.core import FlowsheetBlock
from proteuslib.unit_models.reverse_osmosis_0D import ReverseOsmosis0D, ConcentrationPolarizationType, MassTransferCoefficient, PressureChangeType
# from proteuslib.unit_models.reverse_osmosis_1D import ReverseOsmosis1D, PressureChangeType, ConcentrationPolarizationType, MassTransferCoefficient

# Import log10 function from pyomo
from pyomo.environ import log10

import idaes.logger as idaeslog

def dummy_method(b, *args, **kwargs):
    return 42

# TODO:

def rule_diffus_phase_comp(b, p, j):  # diffusivity, eq 6 in Bartholomew
    # diffus_param_dict = {'0': 1.51e-9, '1': -2.00e-9, '2': 3.01e-8,
    #                      '3': -1.22e-7, '4': 1.53e-7}
    #
    # diffus_param = Param(
    #     diffus_param_dict.keys(),
    #     domain=Reals,
    #     initialize=diffus_param_dict,
    #     units=pyunits.m ** 2 / pyunits.s,
    #     doc='Dynamic viscosity parameters')
    mw = b.params.get_component(j).mw
    mass_frac_phase_comp = b.conc_mol_phase_comp['Liq', j] * mw / b.dens_mass_phase['Liq']
    return b.diffus_phase_comp[p, j] == (1.53e-7* mass_frac_phase_comp ** 4
                                         + (-1.22e-7) * mass_frac_phase_comp ** 3
                                         + (3.01e-8) * mass_frac_phase_comp ** 2
                                         + (-1.22e-7) * mass_frac_phase_comp
                                         + 1.53e-7)

# Configuration dictionary
water_thermo_config = {
    "components": {
        'H2O': {"type": Solvent,
              # Define the methods used to calculate the following properties
              "dens_mol_liq_comp": Perrys,
              # "enth_mol_liq_comp": Perrys,
              "cp_mol_liq_comp": Perrys,
              "entr_mol_liq_comp": Perrys,
              "relative_permittivity_liq_comp": relative_permittivity_constant,
              # Parameter data is always associated with the methods defined above
              "parameter_data": {
                    "mw": (18.0153, pyunits.g/pyunits.mol),
                    "relative_permittivity_liq_comp": 78.54,
                    "pressure_crit": (220.64E5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    # Comes from Perry's Handbook:  p. 2-98
                    "dens_mol_liq_comp_coeff": {
                        '1': (5.459, pyunits.kmol*pyunits.m**-3),
                        '2': (0.30542, pyunits.dimensionless),
                        '3': (647.13, pyunits.K),
                        '4': (0.081, pyunits.dimensionless)},
                    "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ/pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (0, pyunits.kJ/pyunits.mol),
                    # Comes from Perry's Handbook:  p. 2-174
                    "cp_mol_liq_comp_coeff": {
                        '1': (2.7637E5, pyunits.J/pyunits.kmol/pyunits.K),
                        '2': (-2.0901E3, pyunits.J/pyunits.kmol/pyunits.K**2),
                        '3': (8.125, pyunits.J/pyunits.kmol/pyunits.K**3),
                        '4': (-1.4116E-2, pyunits.J/pyunits.kmol/pyunits.K**4),
                        '5': (9.3701E-6, pyunits.J/pyunits.kmol/pyunits.K**5)},
                    "cp_mol_ig_comp_coeff": {
                        'A': (30.09200, pyunits.J/pyunits.mol/pyunits.K),
                        'B': (6.832514, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-1),
                        'C': (6.793435, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-2),
                        'D': (-2.534480, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**-3),
                        'E': (0.082139, pyunits.J*pyunits.mol**-1*pyunits.K**-1*pyunits.kiloK**2),
                        'F': (-250.8810, pyunits.kJ/pyunits.mol),
                        'G': (223.3967, pyunits.J/pyunits.mol/pyunits.K),
                        'H': (0, pyunits.kJ/pyunits.mol)},
                    "entr_mol_form_liq_comp_ref": (69.95, pyunits.J/pyunits.K/pyunits.mol),
                    "pressure_sat_comp_coeff": {
                        'A': (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                        'B': (1435.264, pyunits.K),
                        'C': (-64.848, pyunits.K)}
                                },
                    # End parameter_data
                    },
        "TDS": {"type": Solute,
                 # "enth_mol_liq_comp": dummy_method,
                 "dens_mol_liq_comp": dummy_method,
                 "parameter_data": {
                                    "mw": 58.44e-3 * pyunits.kg/ pyunits.mol}
                 },
        # "CaSO4": {"type": Apparent,
        #          "dissociation_species": {"Ca_2+": 1, "SO4_2-": 1}},
        # "Ca_2+": {"type": Solute,
        #           },
        # "SO4_2-": {"type": Solute,
        #           },
        "Na+": {"type": Cation,
                "charge": +1,
                "parameter_data": {
                    "mw": 23e-3 * pyunits.kg / pyunits.mol}
                # "enth_mol_liq_comp": dummy_method
                },
        "Cl-": {"type": Anion,
                "charge": -1,
                "parameter_data": {
                    "mw": 35e-3 * pyunits.kg / pyunits.mol}
                # "enth_mol_liq_comp": dummy_method
               },
        # "Mg_2+": {"type": Solute,
        #         }
                },  # End Component list
        "phases":  {'Liq': {"type": AqueousPhase,
                            "equation_of_state": ENRTL,
                            "equation_of_state_options": {
                                    "reference_state": Unsymmetric},
                            "diffus_phase_comp": rule_diffus_phase_comp
                            "visc_d_phase":

                            },
                    },
        "state_definition": FTPx,
        "state_bounds": {"flow_mol": (0, 50, 100),
                         "temperature": (273.15, 300, 650),
                         "pressure": (5e4, 1e5, 1e6)
                     },
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},
        "state_components": StateIndex.true,
}
    # End thermo_config definition

m = ConcreteModel()
m.fs = FlowsheetBlock(default={'dynamic': False})
m.fs.properties = GenericParameterBlock(default=water_thermo_config)
m.fs.ro = ReverseOsmosis0D(default={"property_package": m.fs.properties,
                                    "has_pressure_change": True,
                                    "concentration_polarization_type": ConcentrationPolarizationType.calculated,
                                    "mass_transfer_coefficient": MassTransferCoefficient.calculated,
                                    "pressure_change_type": PressureChangeType.calculated})
