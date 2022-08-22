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

import numpy as np
import sys, os
import matplotlib.pyplot as plt

from pyomo.environ import (
    Param,
    Expression,
    Objective,
    SolverFactory,
    units as pyunits,
    Var,
)

from pyomo.util.check_units import assert_units_consistent

from watertap.examples.flowsheets.RO_multiperiod_model.multiperiod_RO import (
    create_multiperiod_swro_model,
)

from idaes.core.solvers import get_solver


def main(
    ndays=1,
    filename="dagget_CA_LMP_hourly_2015.csv",
):
    file_path = os.path.realpath(__file__)
    base_path = os.path.dirname(file_path)
    data_path = os.path.join(base_path, filename)
    # number of time steps assuming 1 hour as the base time
    n_steps = int(ndays * 24)

    # get data
    data = _get_lmp(n_steps, data_path)

    mp_swro = build_flowsheet(n_steps)

    m, t_blocks = set_objective(mp_swro, data)

    # TODO: uncomment after checking the results for lack of variability in flowrate/pressure
    # , then resolve units error
    # assert_units_consistent(m)
    m, _ = solve(m)

    return m, t_blocks, data


def _get_lmp(time_steps, data_path):
    """
    Get price signals from data set
    :param time_steps: Number of time steps considered in MP analysis
    :param data_path: Price [$/kWh] on the same interval as time_steps
    :return: reshaped data
    """
    # read data into array from file
    lmp_data = np.genfromtxt(data_path, delimiter=",")

    # index only the desired number of timesteps
    return lmp_data[:time_steps]


def build_flowsheet(n_steps):
    # create mp model
    mp_swro = create_multiperiod_swro_model(n_time_points=n_steps)

    return mp_swro


def set_objective(mp_swro, lmp, carbontax=3.5):
    # Retrieve pyomo model and active process blocks (i.e. time blocks)
    m = mp_swro.pyomo_model
    t_blocks = mp_swro.get_active_process_blocks()

    # index the flowsheet for each timestep
    for count, blk in enumerate(t_blocks):
        blk_swro = blk.ro_mp

        # set price and carbon signals as parameters
        blk.lmp_signal = Param(
            default=lmp[count],
            mutable=True,
            units=blk_swro.fs.costing.base_currency / pyunits.kWh,
        )
        blk.carbon_intensity = Param(
            default=100, mutable=True, units=pyunits.kg / pyunits.MWh
        )
        blk.carbon_tax = Param(
            default=carbontax,
            mutable=True,
            units=blk_swro.fs.costing.base_currency / pyunits.kg,
        )

        # set the electricity_price in each flowsheet
        blk_swro.fs.costing.electricity_cost.fix(blk.lmp_signal)

        # combine/place flowsheet level cost metrics on each time block
        blk.weighted_LCOW = Expression(
            expr=blk_swro.fs.costing.LCOW * blk_swro.fs.costing.annual_water_production,
            doc="annual flow weighted LCOW",
        )
        blk.water_prod = Expression(
            expr=blk_swro.fs.costing.annual_water_production,
            doc="annual water production",
        )
        blk.energy_consumption = Expression(
            expr=blk_swro.fs.costing.specific_energy_consumption
            * pyunits.convert(
                blk_swro.fs.costing.annual_water_production,
                to_units=pyunits.m**3 / pyunits.hour,
            ),
            doc="Energy consumption per timestep ",
        )
        blk.carbon_emission = Expression(
            expr=blk.energy_consumption
            * pyunits.convert(blk.carbon_intensity, to_units=pyunits.kg / pyunits.kWh),
            doc="Equivalent carbon emissions per timestep ",
        )
        blk.annual_carbon_cost = Expression(
            expr=blk.carbon_tax
            * pyunits.convert(blk.carbon_emission, to_units=pyunits.kg / pyunits.year)
            * (blk_swro.fs.costing.base_currency / pyunits.kg),
            doc="Annual cost associated with carbon emissions, to be used in the objective",
        )

        # deactivate each block-level objective function
        blk_swro.fs.objective.deactivate()

    # compile time block level expressions into a model-level objective
    m.obj = Objective(
        expr=(
            sum([blk.weighted_LCOW for blk in t_blocks])
            + sum([blk.annual_carbon_cost for blk in t_blocks])
        )
        / sum([blk.water_prod for blk in t_blocks]),
        doc="Flow-integrated average cost and carbon tax on an annual basis",
    )

    # fix the initial pressure to default operating pressure at 1 kg/s and 50% recovery
    t_blocks[0].ro_mp.previous_pressure.fix(55e5)

    return m, t_blocks


def solve(m):
    # solve
    # opt = SolverFactory("ipopt")
    opt = get_solver()
    results = opt.solve(m, tee=True)
    print("Solver:", opt.__class__)
    return m, results


def visualize_results(m, t_blocks, data):
    time_step = np.array(range(len(t_blocks)))
    recovery = np.array(
        [
            blk.ro_mp.fs.RO.recovery_mass_phase_comp[0, "Liq", "H2O"].value
            for blk in t_blocks
        ]
    )

    pump1_flow = np.array(
        [
            blk.ro_mp.fs.P1.control_volume.properties_out[0].flow_vol_phase["Liq"]()
            for blk in t_blocks
        ]
    )
    pump2_flow = np.array(
        [
            blk.ro_mp.fs.P2.control_volume.properties_out[0].flow_vol_phase["Liq"]()
            for blk in t_blocks
        ]
    )
    pressure = np.array(
        [
            blk.ro_mp.fs.P1.control_volume.properties_out[0].pressure.value
            for blk in t_blocks
        ]
    )
    power = np.array([blk.energy_consumption() for blk in t_blocks])
    pump1_efficiency = np.array(
        [blk.ro_mp.fs.P1.efficiency_pump[0]() for blk in t_blocks]
    )
    pump2_efficiency = np.array(
        [blk.ro_mp.fs.P2.efficiency_pump[0]() for blk in t_blocks]
    )

    fig, ax = plt.subplots(3, 2)
    fig.suptitle(
        "SWRO MultiPeriod optimization results: {} hours\nLCOW: {} $/m3".format(
            len(t_blocks), round(m.obj(), 2)
        )
    )

    ax[0, 0].plot(time_step, data)
    ax[0, 0].set_xlabel("Time [hr]")
    ax[0, 0].set_ylabel("Electricity price [$/kWh]")

    ax[1, 0].plot(time_step, pump1_flow * 3600, label="P1")
    ax[1, 0].plot(time_step, pump2_flow * 3600, label="P2")
    ax[1, 0].set_xlabel("Time [hr]")
    ax[1, 0].set_ylabel("Pump flowrate [m3/hr]")

    ax[2, 0].plot(time_step, recovery * 100)
    ax[2, 0].set_xlabel("Time [hr]")
    ax[2, 0].set_ylabel("Water Recovery [%]")

    ax[0, 1].plot(time_step, power)
    ax[0, 1].set_xlabel("Time [hr]")
    ax[0, 1].set_ylabel("Net Power [kWh]")

    ax[1, 1].plot(time_step, pump1_efficiency * 100, label="P1")
    ax[1, 1].plot(time_step, pump2_efficiency * 100, label="P2")
    ax[1, 1].set_xlabel("Time [hr]")
    ax[1, 1].set_ylabel("Pump efficiency [%]")
    ax[1, 1].legend()

    ax[2, 1].plot(time_step, pressure * 1e-5)
    ax[2, 1].set_xlabel("Time [hr]")
    ax[2, 1].set_ylabel("RO Inlet Pressure [bar]")


if __name__ == "__main__":
    m, t_blocks, data = main(*sys.argv[1:])
    visualize_results(m, t_blocks, data)
