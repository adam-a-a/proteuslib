#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 12:14:49 2022

@author: Maggie
"""

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Expression,
    value,
    Block,
    Constraint,
    Var,
    Param,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import TransformationFactory, assert_optimal_termination
from pyomo.network import Arc

# Imports from IDAES
# Import flowsheet block from IDAES core
from idaes.core import FlowsheetBlock

# Import function to get default solver
from idaes.core.solvers import get_solver
from idaes.core.util.initialization import propagate_state

# Import function to check degrees of freedom
from idaes.core.util.model_statistics import degrees_of_freedom

# Import utility function for calculating scaling factors
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import Product
from idaes.core import UnitModelCostingBlock
from watertap.unit_models.zero_order import IonExchangeReactiveZO
from watertap.unit_models.zero_order import NanofiltrationZO
from watertap.unit_models.zero_order import PumpElectricityZO  # PumpZO not in costing
from watertap.unit_models.zero_order import ChemicalAdditionZO
from watertap.unit_models.zero_order import FeedZO

from watertap.core.wt_database import Database
from watertap.core.zero_order_properties import WaterParameterBlock
from watertap.core.zero_order_costing import ZeroOrderCosting
from watertap.core.util import assert_degrees_of_freedom

solver = get_solver()


def main():
    m = build()

    set_operating_conditions(m)
    assert_degrees_of_freedom(m, 0)
    assert_units_consistent(m)

    initialize_system(m)

    results = solve(m)
    assert_optimal_termination(results)
    # display_results(m.fs)
    add_costing(m)
    # initialize costing
    m.fs.costing.initialize()

    assert_degrees_of_freedom(m, 0)
    results = solve(m)
    assert_optimal_termination(results)
    # display_costing(m)
    return m, results


def build():
    # TODO: Over-specified rn and getting error of "No variables appear in
    # the Pyomo model constraints or objective.
    # This is not supported by the NL file interface" even when not over-specified
    # setup model
    m = ConcreteModel()
    m.db = Database()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = WaterParameterBlock(default={"solute_list": ["tds", "calcium"]})
    m.fs.costing = ZeroOrderCosting()

    # create unit models
    m.fs.feed = FeedZO(default={"property_package": m.fs.properties})
    m.fs.pump1 = PumpElectricityZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.pump2 = PumpElectricityZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.pump3 = PumpElectricityZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.ion_ex = IonExchangeReactiveZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.nf1 = NanofiltrationZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.nf2 = NanofiltrationZO(
        default={"property_package": m.fs.properties, "database": m.db}
    )
    m.fs.nacl_addition = ChemicalAdditionZO(
        default={
            "property_package": m.fs.properties,
            "database": m.db,
            "process_subtype": "sodium_chloride",  # ADD THIS TO YAML! (did)
        }
    )
    m.fs.product = Product(default={"property_package": m.fs.properties})
    # m.fs.disposal = Product(default={"property_package": m.fs.properties})

    # connect unit models
    m.fs.s1 = Arc(source=m.fs.feed.outlet, destination=m.fs.nacl_addition.inlet)
    m.fs.s2 = Arc(source=m.fs.nacl_addition.outlet, destination=m.fs.pump1.inlet)
    m.fs.s3 = Arc(source=m.fs.pump1.outlet, destination=m.fs.ion_ex.inlet)
    m.fs.s4 = Arc(source=m.fs.ion_ex.byproduct, destination=m.fs.pump2.inlet)
    m.fs.s5 = Arc(source=m.fs.pump2.outlet, destination=m.fs.nf1.inlet)
    m.fs.s6 = Arc(source=m.fs.nf1.treated, destination=m.fs.pump3.inlet)
    m.fs.s7 = Arc(source=m.fs.pump3.outlet, destination=m.fs.nf2.inlet)
    m.fs.s8 = Arc(source=m.fs.nf2.treated, destination=m.fs.product.inlet)
    # m.fs.s9 = Arc(source=m.fs.ion_ex.treated, destination=m.fs.disposal.inlet)
    # m.fs.s10 = Arc(source=m.fs.nf1.byproduct, destination=m.fs.disposal.inlet)
    # m.fs.s11 = Arc(source=m.fs.nf2.byproduct, destination=m.fs.disposal.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m)

    # scaling (NEED TO FIGURE THIS OUT?? SEEMS LIKE SHOULD HAVE MORE)
    iscale.calculate_scaling_factors(m.fs)

    return m


# specify the flowsheet
def set_operating_conditions(m):
    # m.fs.feed.flow_vol[0].fix(1)
    # m.fs.feed.conc_mass_comp[0, "tds"].fix(60)
    # m.fs.feed.conc_mass_comp[0, "calcium"].fix(1)

    m.fs.feed.properties[0].flow_mass_comp["H2O"].fix(265)
    m.fs.feed.properties[0].flow_mass_comp["tds"].fix(0.6625)
    m.fs.feed.properties[0].flow_mass_comp["calcium"].fix(0.219)
    solve(m.fs.feed)

    # m.db.get_unit_operation_parameters("pump_electricity")
    # m.fs.pump1.load_parameters_from_database()
    # m.fs.pump2.load_parameters_from_database()
    # m.fs.pump3.load_parameters_from_database()
    m.fs.pump1.applied_pressure.fix(1.01325)
    m.fs.pump1.eta_pump.fix(0.8)
    m.fs.pump1.eta_motor.fix(0.8)
    m.fs.pump2.applied_pressure.fix(5)
    m.fs.pump2.eta_pump.fix(0.8)
    m.fs.pump2.eta_motor.fix(0.8)
    m.fs.pump3.applied_pressure.fix(5)
    m.fs.pump3.eta_pump.fix(0.8)
    m.fs.pump3.eta_motor.fix(0.8)

    # m.db.get_unit_operation_parameters("nanofiltration")
    m.fs.nf1.load_parameters_from_database(use_default_removal=True)
    m.fs.nf2.load_parameters_from_database(use_default_removal=True)
    m.fs.nf1.recovery_frac_mass_H2O.fix(0.8)
    m.fs.nf2.recovery_frac_mass_H2O.fix(0.8)

    m.db.get_unit_operation_parameters("ion_exchange")
    m.fs.ion_ex.load_parameters_from_database(use_default_removal=True)

    # m.db.get_unit_operation_parameters("chemical_addition")
    # m.fs.nacl_addition.load_parameters_from_database()
    m.fs.nacl_addition.chemical_dosage.fix(10)
    m.fs.nacl_addition.solution_density.fix(9.57)
    m.fs.nacl_addition.ratio_in_solution.fix(1)

    print("DOF = ", degrees_of_freedom(m))
    assert_units_consistent(m)


def initialize_system(m):
    # solve feed
    # solver.solve(m.fs.feed) (gave an error)
    m.fs.feed.initialize()

    # initialize units
    propagate_state(m.fs.s1)
    m.fs.nacl_addition.initialize()

    propagate_state(m.fs.s2)
    m.fs.pump1.initialize()

    propagate_state(m.fs.s3)
    m.fs.ion_ex.initialize()

    propagate_state(m.fs.s4)
    m.fs.pump2.initialize()

    propagate_state(m.fs.s5)
    m.fs.nf1.initialize()

    propagate_state(m.fs.s6)
    m.fs.pump3.initialize()

    propagate_state(m.fs.s7)
    m.fs.nf2.initialize()

    propagate_state(m.fs.s8)
    # propagate_state(m.fs.s9)
    # propagate_state(m.fs.s10)
    # propagate_state(m.fs.s11)


def solve(m):
    # solve model
    results = solver.solve(m, tee=True)
    # print(results)
    return results


def add_costing(m):
    # unit equipment capital and operating costs
    m.fs.pump1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.pump2.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.pump3.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.ion_ex.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.nf1.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.nf2.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )
    m.fs.nacl_addition.costing = UnitModelCostingBlock(
        default={"flowsheet_costing_block": m.fs.costing}
    )

    # system costing - total investment and operating costs
    m.fs.costing.cost_process()
    m.fs.costing.add_electricity_intensity(m.fs.product.properties[0].flow_vol)
    m.fs.costing.add_LCOW(m.fs.product.properties[0].flow_vol)


if __name__ == "__main__":
    m, res = main()
