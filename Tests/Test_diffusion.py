import unittest
import TB_Model

def parameter_setup():
    params = dict()
    params['max_depth'] = 3
    params['initial_oxygen'] = 1
    params['blood_vessel_value'] = 1.5
    params['time_step'] = 0.001
    params['spatial_step'] = 0.2
    params['oxygen_diffusion'] = 1
    params['oxygen_from_source'] = 2.4
    params['oxygen_uptake_from_bacteria'] = 1
    params['chemotherapy_diffusion'] = 0.75
    params['chemotherapy_from_source'] = 1
    params['chemotherapy_decay'] = 0.35
    params['chemokine_diffusion'] = 0.05
    params['chemokine_from_bacteria'] = 0
    params['chemokine_from_macrophage'] = 1
    params['chemokine_decay'] = 0.347
    params['caseum_distance_to_reduce_diffusion'] = 2
    params['caseum_threshold_to_reduce_diffusion'] = 2
    params['oxygen_diffusion_caseum_reduction'] = 1.5
    params['chemotherapy_diffusion_caseum_reduction'] = 1.5
    params['chemotherapy_schedule1_start_lower'] = 168
    params['chemotherapy_schedule1_start_upper'] = 336
    params['chemotherapy_schedule1_end'] = 700
    params['chemotherapy_schedule2_start'] = 900
    params['bacteria_replication_fast_upper'] = 31
    params['bacteria_replication_fast_lower'] = 17
    params['bacteria_replication_slow_upper'] = 72
    params['bacteria_replication_slow_lower'] = 48
    params['bacteria_threshold_for_t_cells'] = 50
    params['t_cell_recruitment_probability'] = 2
    params['chemokine_scale_for_t_cell_recruitment'] = 0.01
    params['bacteria_threshold_for_macrophage_recruitment'] = 70
    params['macrophage_recruitment_probability'] = 7
    params['chemokine_scale_for_macrophage_recruitment_above_threshold'] = 0.001
    params['chemokine_scale_for_macrophage_recruitment_below_threshold'] = 0.1
    params['chemotherapy_scale_for_kill_fast_bacteria'] = 10
    params['chemotherapy_scale_for_kill_slow_bacteria'] = 20
    params['chemotherapy_scale_for_kill_macrophage'] = 30
    params['t_cell_movement_time'] = 167
    params['t_cell_age_threshold'] = 72000
    params['t_cell_random_move_probability'] = 20
    params['t_cell_kills_macrophage_probability'] = 75
    params['resting_macrophage_age_limit'] = 2400000
    params['resting_macrophage_movement_time'] = 333
    params['prob_resting_macrophage_random_move'] = 10
    params['minimum_chemokine_for_resting_macrophage_movement'] = 20
    params['active_macrophage_age_limit'] = 240000
    params['active_macrophage_movement_time'] = 7800
    params['prob_active_macrophage_kill_fast_bacteria'] = 20
    params['prob_active_macrophage_kill_slow_bacteria'] = 30
    params['infected_macrophage_age_limit'] = 2400000
    params['infected_macrophage_movement_time'] = 240000
    params['chronically_infected_macrophage_age_limit'] = 2400000
    params['chronically_infected_macrophage_movement_time'] = 240000
    params['chemokine_scale_for_macrophage_activation'] = 50
    params['chemokine_scale_for_macrophage_deactivation'] = 1
    params['bacteria_to_turn_chronically_infected'] = 10
    params['bacteria_to_burst_macrophage'] = 20
    params['oxygen_scale_for_metabolism_change_to_slow'] = 1
    params['oxygen_scale_for_metabolism_change_to_fast'] = 99
    params['interval_to_record_results'] = 1000

    return params

class OxygenDiffusionTestCase(unittest.TestCase):

    def setUp(self):
        params = parameter_setup()

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.params = params
        self.atts = atts

        blood_vessels = [(25, 25)]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [100, 100], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)
        self.sort_out_halo()

    def tearDown(self):
        for a in self.topology.automata:
            a.close_files()

    def sort_out_halo(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_something(self):
        for a in self.topology.automata:
            a.diffusion_pre_process()

        for a in self.topology.automata:
            a.diffusion(False)
            a.swap_grids()

        self.assertAlmostEqual(self.topology.automata[0].grid[(25, 25)]['oxygen'],  1.3536)
        self.assertAlmostEqual(self.topology.automata[0].grid[(24, 25)]['oxygen'],  0.0375)
        self.assertAlmostEqual(self.topology.automata[0].grid[(25, 24)]['oxygen'],  0.0375)
        self.assertAlmostEqual(self.topology.automata[0].grid[(25, 26)]['oxygen'],  0.0375)
        self.assertAlmostEqual(self.topology.automata[0].grid[(26, 25)]['oxygen'],  0.0375)



if __name__ == '__main__':
    unittest.main()
