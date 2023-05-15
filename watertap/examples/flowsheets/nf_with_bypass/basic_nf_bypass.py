#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

from pyomo.environ import (
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core import (
    FlowsheetBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.models.unit_models import (
    Mixer,
    Separator,
    Product,
    Feed,
    StateJunction,
)
from idaes.models.unit_models.mixer import MomentumMixingType, MixingType

from watertap.unit_models.nanofiltration_DSPMDE_0D import (
    NanofiltrationDSPMDE0D,
)

from watertap.unit_models.pressure_changer import Pump

from pyomo.environ import ConcreteModel, TransformationFactory

import math

from watertap.property_models.multicomp_aq_sol_prop_pack import (
    MCASParameterBlock,
    ActivityCoefficientModel,
    DensityCalculation,
)
import idaes.core.util.scaling as iscale


def main():
    m = build()
    feed_mass_frac = {
        "Ca_2+": 4.0034374454637006e-03,
        "HCO3_-": 0.00022696833343821863,
        "SO4_2-": 0.00020497140244420624,
        "Cl_-": 0.0004559124032433401,
        "Na_+": 0.00043333830389924205,
    }

    mole_frac = {
        "Ca_2+": 0.0007992007992007992,
        "HCO3_-": 0.0029762294974498824,
        "SO4_2-": 0.0017072662919008952,
        "Cl_-": 0.010289138002902268,
        "Na_+": 0.015081448209232778,
        "H2O": 44.39013800000001,
    }

    set_NF_feed(m.fs, feed_mass_frac=None, mole_frac=mole_frac)
    init(m)

    solver = get_solver()
    solver.solve(m, tee=True)


def build():

    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    default = {
        "solute_list": [
            "Ca_2+",
            "SO4_2-",
            "HCO3_-",
            "Na_+",
            "Cl_-",
        ],
        "diffusivity_data": {
            ("Liq", "Ca_2+"): 9.2e-10,
            ("Liq", "SO4_2-"): 1.06e-9,
            ("Liq", "HCO3_-"): 1.19e-9,
            ("Liq", "Na_+"): 1.33e-9,
            ("Liq", "Cl_-"): 2.03e-9,
        },
        "mw_data": {
            "H2O": 18e-3,
            "Ca_2+": 40e-3,
            "HCO3_-": 61.0168e-3,
            "SO4_2-": 96e-3,
            "Na_+": 23e-3,
            "Cl_-": 35e-3,
        },
        "stokes_radius_data": {
            "Ca_2+": 0.309e-9,
            "HCO3_-": 2.06e-10,
            "SO4_2-": 0.230e-9,
            "Cl_-": 0.121e-9,
            "Na_+": 0.184e-9,
        },
        "charge": {
            "Ca_2+": 2,
            "HCO3_-": -1,
            "SO4_2-": -2,
            "Na_+": 1,
            "Cl_-": -1,
        },
        "activity_coefficient_model": ActivityCoefficientModel.ideal,
        "density_calculation": DensityCalculation.constant,
    }

    m.fs.properties = MCASParameterBlock(**default)

    m.fs.feed = Feed(property_package=m.fs.properties)
    m.fs.product = Product(property_package=m.fs.properties)
    m.fs.by_pass_splitter = Separator(
        property_package=m.fs.properties,
        outlet_list=["nf_stage", "bypass"],
    )

    m.fs.NF = FlowsheetBlock(dynamic=False)
    m.fs.NF.feed = StateJunction(property_package=m.fs.properties)
    m.fs.NF.product = StateJunction(property_package=m.fs.properties)

    m.fs.NF.pump = Pump(property_package=m.fs.properties)
    m.fs.NF.nfUnit = NanofiltrationDSPMDE0D(property_package=m.fs.properties)

    m.fs.NF.feed_to_pump = Arc(
        source=m.fs.NF.feed.outlet, destination=m.fs.NF.pump.inlet
    )
    m.fs.NF.pump_to_nf = Arc(
        source=m.fs.NF.pump.outlet, destination=m.fs.NF.nfUnit.inlet
    )
    m.fs.NF.nf_to_product = Arc(
        source=m.fs.NF.nfUnit.permeate, destination=m.fs.NF.product.inlet
    )

    m.fs.total_product_mixer = Mixer(
        property_package=m.fs.properties,
        inlet_list=["bypass", "nf_stage"],
        energy_mixing_type=MixingType.none,
        momentum_mixing_type=MomentumMixingType.minimize,
    )

    m.fs.feed_to_splitter = Arc(
        source=m.fs.feed.outlet, destination=m.fs.by_pass_splitter.inlet
    )

    m.fs.splitter_to_nffeed = Arc(
        source=m.fs.by_pass_splitter.nf_stage, destination=m.fs.NF.feed.inlet
    )

    m.fs.splitter_to_mixer = Arc(
        source=m.fs.by_pass_splitter.bypass, destination=m.fs.total_product_mixer.bypass
    )

    m.fs.nfUnit_product_to_mixer = Arc(
        source=m.fs.NF.product.outlet,
        destination=m.fs.total_product_mixer.nf_stage,
    )
    m.fs.mixer_to_product = Arc(
        source=m.fs.total_product_mixer.outlet, destination=m.fs.product.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)
    return m


def init(m):
    solver = get_solver()
    m.fs.feed.initialize(optarg=solver.options)
    m.fs.feed.properties[0].pressure.fix(101325)
    propagate_state(m.fs.feed_to_splitter)

    m.fs.by_pass_splitter.split_fraction[0, "bypass"].fix(0.5)
    m.fs.by_pass_splitter.mixed_state.initialize(optarg=solver.options)
    m.fs.by_pass_splitter.initialize(optarg=solver.options)

    propagate_state(m.fs.splitter_to_nffeed)
    m.fs.NF.feed.initialize(optarg=solver.options)
    propagate_state(m.fs.NF.feed_to_pump)
    m.fs.NF.pump.outlet.pressure[0].fix(2.0 * 1e5)
    m.fs.NF.pump.efficiency_pump[0].fix(0.75)

    iscale.set_scaling_factor(m.fs.NF.pump.control_volume.work, 1e-4)

    m.fs.NF.pump.initialize(optarg=solver.options)

    propagate_state(m.fs.NF.pump_to_nf)

    m.fs.NF.nfUnit.recovery_vol_phase.fix(0.1)
    m.fs.NF.nfUnit.area.unfix()
    m.fs.NF.nfUnit.area = 50

    m.fs.NF.nfUnit.spacer_porosity.fix(0.85)
    m.fs.NF.nfUnit.channel_height.fix(1e-3)
    m.fs.NF.nfUnit.velocity[0, 0].fix(0.25)
    m.fs.NF.nfUnit.spacer_mixing_efficiency.fix()
    m.fs.NF.nfUnit.spacer_mixing_length.fix()

    m.fs.NF.nfUnit.radius_pore.fix(0.5e-9)

    m.fs.NF.nfUnit.membrane_thickness_effective.fix(8.598945196055952e-07)
    m.fs.NF.nfUnit.membrane_charge_density.fix(-680)
    m.fs.NF.nfUnit.dielectric_constant_pore.fix(41.3)
    m.fs.NF.nfUnit.mixed_permeate[0].pressure.fix(101325)

    m.fs.NF.nfUnit.initialize(
        optarg=solver.options,
        automate_rescale=False,
        initialize_guess={
            "solvent_recovery": 0.01,
            "solute_recovery": 0.01,
            "cp_modulus": 1,
        },
    )
    propagate_state(m.fs.NF.nf_to_product)
    m.fs.NF.product.initialize(optarg=solver.options)
    propagate_state(m.fs.nfUnit_product_to_mixer)
    propagate_state(m.fs.splitter_to_mixer)
    m.fs.total_product_mixer.mixed_state.initialize(optarg=solver.options)
    m.fs.total_product_mixer.initialize(optarg=solver.options)


def set_NF_feed(blk, feed_mass_frac=None, mole_frac=None):
    mass_flow_in = 1 * pyunits.kg / pyunits.s
    if feed_mass_frac is not None:
        for ion, x in feed_mass_frac.items():
            mol_comp_flow = (
                x
                * pyunits.kg
                / pyunits.kg
                * mass_flow_in
                / blk.feed.properties[0].mw_comp[ion]
            )
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(mol_comp_flow)
        H2O_mass_frac = 1 - sum(x for x in feed_mass_frac.values())
        H2O_mol_comp_flow = (
            H2O_mass_frac
            * pyunits.kg
            / pyunits.kg
            * mass_flow_in
            / blk.feed.properties[0].mw_comp["H2O"]
        )
        blk.feed.properties[0].flow_mol_phase_comp["Liq", "H2O"].fix(H2O_mol_comp_flow)

    if mole_frac is not None:
        for ion, x in mole_frac.items():
            blk.feed.properties[0].flow_mol_phase_comp["Liq", ion].fix(x)

    set_NF_feed_scaling(blk)
    blk.feed.properties[0].assert_electroneutrality(
        defined_state=True,
        adjust_by_ion="Cl_-",
        get_property="flow_mol_phase_comp",
    )

    blk.feed.properties[0].temperature.fix(298.15)
    iscale.calculate_scaling_factors(blk)


def calc_scale(value):
    return -1 * math.floor(math.log(value, 10))


def set_NF_feed_scaling(blk):
    _add = 0
    for index in blk.feed.properties[0].flow_mol_phase_comp:
        scale = calc_scale(blk.feed.properties[0].flow_mol_phase_comp[index].value)
        print(f"{index} flow_mol_phase_comp scaling factor = {10**(scale+_add)}")
        blk.properties.set_default_scaling(
            "flow_mol_phase_comp", 10 ** (scale + 0), index=index
        )


if __name__ == "__main__":
    main()