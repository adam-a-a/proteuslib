
import os
import proteuslib.flowsheets.full_treatment_train.analysis.flowsheet_single_stage as fs_single_stage
from proteuslib.tools.parameter_sweep import _init_mpi, LinearSample, parameter_sweep
import proteuslib.flowsheets.full_treatment_train.example_flowsheets.costing as costing

def append_costing_outputs(m, outputs=None):
    costing.display_cost_breakdown(m)
    for unit_name in units_to_cost:
        for cost_type in ['capital_cost', 'operating_cost']:
            unit = getattr(m.fs, unit_name)
            cost = getattr(unit.costing, cost_type)

            outputs['%s_%s' % (cost_type, unit_name)] = cost

    return outputs

def run_analysis(nx, RO_type):
    desal_kwargs = {'has_desal_feed': False, 'is_twostage': False, 'has_ERD': True,
                    'RO_type': RO_type, 'RO_base': 'TDS', 'RO_level': 'detailed'}

    m = fs_single_stage.optimize_flowsheet(**desal_kwargs)


    sweep_params = {}
    sweep_params['System Recovery'] = LinearSample(m.fs.system_recovery_target, 0.3, 0.6, nx)

    outputs = {}
    outputs['LCOW'] = m.fs.costing.LCOW
    outputs['Saturation Index'] = m.fs.desal_saturation.saturation_index
    outputs['Pump Pressure'] = m.fs.pump_RO.control_volume.properties_out[0].pressure
    outputs, _, _ = costing.display_cost_breakdown(m, cost_type='levelized')
    output_filename = 'output/fs_single_stage/results_%sRO.csv' % (desal_kwargs['RO_type'])

    opt_function = fs_single_stage.optimize

    global_results = parameter_sweep(m, sweep_params, outputs, output_filename,
                                     optimize_function=opt_function,
                                     debugging_data_dir=os.path.split(output_filename)[0]+'/local')
    return global_results, sweep_params

if __name__ == "__main__":
    nx= 4
    RO_type= '0D'
    global_results, sweep_params = run_analysis(nx, RO_type)
    print(global_results)