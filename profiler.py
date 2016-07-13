import cProfile
import TB_Model


def setup_topology():

    params = dict()
    params['max_depth'] = 3
    params['caseum_threshold_to_reduce_diffusion'] = 999
    params['initial_oxygen'] = 1.5
    params['oxygen_diffusion'] = 0.0
    params['chemotherapy_diffusion'] = 0.0
    params['caseum_distance_to_reduce_diffusion'] = 2
    params['spatial_step'] = 0.2
    params['chemotherapy_schedule1_start'] = 99
    params['chemotherapy_schedule2_start'] = 200
    params['oxygen_from_source'] = 0.0
    params['chemokine_diffusion'] = 0.0
    params['chemokine_decay'] = 0.0
    params['chemokine_from_macrophage'] = 0
    params['bacteria_threshold_for_t_cells'] = 100
    params['chemotherapy_scale_for_kill_macrophage'] = 101
    params['oxygen_uptake_from_bacteria'] = 0
    params['chemokine_from_bacteria'] = 0
    params['bacteria_replication_fast_upper'] = 999
    params['bacteria_replication_fast_lower'] = 998
    params['bacteria_replication_slow_upper'] = 999
    params['bacteria_replication_slow_lower'] = 998
    params['chemotherapy_scale_for_kill_fast_bacteria'] = 101
    params['chemotherapy_scale_for_kill_slow_bacteria'] = 101
    params['resting_macrophage_age_limit'] = 999
    params['active_macrophage_age_limit'] = 999
    params['infected_macrophage_age_limit'] = 999
    params['chronically_infected_macrophage_age_limit'] = 999
    params['resting_macrophage_movement_time'] = 1000
    params['active_macrophage_movement_time'] = 1000
    params['infected_macrophage_movement_time'] = 1000
    params['chronically_infected_macrophage_movement_time'] = 1000
    params['prob_resting_macrophage_random_move'] = 0
    params['minimum_chemokine_for_resting_macrophage_movement'] = 0
    params['bacteria_to_turn_chronically_infected'] = 10
    params['chemokine_scale_for_macrophage_activation'] = 0
    params['chemokine_scale_for_macrophage_deactivation'] = 100
    params['time_step'] = 4
    params['oxygen_scale_for_metabolism_change_to_slow'] = 100
    params['oxygen_scale_for_metabolism_change_to_fast'] = 0
    params['macrophage_recruitment_probability'] = 0

    params['bacteria_to_burst_macrophage'] = 20

    atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
        'chemotherapy', 'chemokine']

    blood_vessels = [[1, 1]]
    fast_bacteria = []
    slow_bacteria = []
    macrophages = [[2, 2]]
    topology = TB_Model.TwoDimensionalTopology([2, 2], [100, 100], atts, params, blood_vessels, fast_bacteria,
                                                    slow_bacteria, macrophages)

    return topology


def sort_out_halos(topology):
    dz = []
    for i in topology.automata:
        dz.append(i.get_danger_zone())
    halos = topology.create_halos(dz)
    for i in range(4):
        topology.automata[i].set_halo(halos[i])

    return topology


def main():
    topology = setup_topology()
    sort_out_halos(topology)

    pr = cProfile.Profile()
    pr.enable()
    topology.automata[0].diffusion_pre_process()
    pr.disable()
    pr.print_stats(sort='cumtime')


if __name__ == '__main__':
    main()