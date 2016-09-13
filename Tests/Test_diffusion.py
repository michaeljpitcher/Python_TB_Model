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


def test_oxygen_diffusion_no_changes(self):
    """
    No diffusion, no oxygen from source. Starts zero, stays zero
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['initial_oxygen'] = 0.0  # init_o2
    self.params['oxygen_from_source'] = 0.0  # gamma[i][j]
    self.params['oxygen_uptake_from_bacteria'] = 0.0  # phi[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt
    self.params['oxygen_diffusion'] = 0.0  # d[i][j]
    self.params['oxygen_diffusion_caseum_reduction'] = 1.0  # value to divide d[i][j] by

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 0.0)


def test_oxygen_diffusion(self):
    """
    Just one blood vessel, no other contents. No value from source so initial oxygen diffuses out, nothing more
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['initial_oxygen'] = 100.0  # init_o2
    self.params['oxygen_from_source'] = 0.0  # gamma[i][j]
    self.params['oxygen_uptake_from_bacteria'] = 0.0  # phi[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt
    self.params['oxygen_diffusion'] = 1.0  # d[i][j]
    self.params['oxygen_diffusion_caseum_reduction'] = 1.0  # value to divide d[i][j] by

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 100.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    # 100 + 0.001 * ((((1+1)/2) * (-100 - 100))/(0.2*0.2) + (((1+1)/2) * (-100 - 100))/(0.2*0.2)) = 90.0
    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 90.0)
    # 0 + 0.001 * ((((1+1)/2) * (100 - 0))/(0.2*0.2) + (((1+1)/2) * (0 - 0))/(0.2*0.2)) = 2.5
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 2.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 2.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 2.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 2.5)

    # Now double the diffusion rate, should result in double the diffusion

    self.params['oxygen_diffusion'] = 2.0  # d[i][j]

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 100.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 80.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 5.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 5.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 5.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 5.0)


def test_oxygen_diffusion_differing_rates(self):
    """
    test oxygen diffusion when the rates differ (due to caseum)
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['initial_oxygen'] = 100.0  # init_o2
    self.params['oxygen_from_source'] = 0.0  # gamma[i][j]
    self.params['oxygen_uptake_from_bacteria'] = 0.0  # phi[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt
    self.params['oxygen_diffusion'] = 1.0  # d[i][j]
    self.params['oxygen_diffusion_caseum_reduction'] = 1.0  # value to divide d[i][j] by

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 100.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    # Manually reduce the diffusion rate
    self.topology.automata[0].grid[2, 1]['oxygen_diffusion_rate'] = 0.5

    self.topology.automata[0].diffusion(False)

    # 100 + 0.001 * (((((1+0.5)/2)*(0 - 100) - ((1+1)/2)*(100 - 0)) / (0.2*0.2)) + ((((1+1)/2)*(0 - 100) - ((1+1)/2)*(100 - 0)) / (0.2*0.2))) = 90.625
    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 90.625)
    # 0 + 0.001 * ((((1+1)/2) * (100 - 0))/(0.2*0.2) + (((1+1)/2) * (0 - 0))/(0.2*0.2)) = 2.5
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 2.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 1.875)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 2.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 2.5)


def test_oxygen_from_source(self):
    """
    No diffusion, only oxygen from source
    :return:
    """
    self.params['blood_vessel_value'] = 1.1  # m[i][j]
    self.params['initial_oxygen'] = 0.0  # init_o2
    self.params['oxygen_from_source'] = 10.0  # gamma[i][j]
    self.params['oxygen_uptake_from_bacteria'] = 0.0  # phi[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt
    self.params['oxygen_diffusion'] = 0.0  # d[i][j]
    self.params['oxygen_diffusion_caseum_reduction'] = 1.0  # value to divide d[i][j] by

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 10 * 1.1 * 0.001)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 0.0)

    # Different oxygen from source rate
    self.params['oxygen_from_source'] = 7.35  # gamma[i][j]

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['oxygen'], 7.35 * 1.1 * 0.001)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['oxygen'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['oxygen'], 0.0)


def test_oxygen_uptake_by_bacteria(self):
    """
    No diffusion, no oxygen from source. oxygen at cell (1,1) is set to a value and decreased by bacteria
    :return:
    """
    self.params['blood_vessel_value'] = 0.0  # m[i][j]
    self.params['initial_oxygen'] = 0.0  # init_o2
    self.params['oxygen_from_source'] = 0.0  # gamma[i][j]
    self.params['oxygen_uptake_from_bacteria'] = 2.0  # phi[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt
    self.params['oxygen_diffusion'] = 0.0  # d[i][j]
    self.params['oxygen_diffusion_caseum_reduction'] = 1.0  # value to divide d[i][j] by

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    bacterium = TB_Model.Bacterium((1, 1), 'fast')
    self.topology.automata[0].bacteria.append(bacterium)
    self.topology.automata[0].grid[1, 1]['contents'] = bacterium
    self.topology.automata[0].grid[1, 1]['oxygen'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 1]['oxygen'], 10.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[1, 1]['oxygen'], 10 - 0.001 * (2 * 10))

    # Change uptake rate

    self.params['oxygen_uptake_from_bacteria'] = 4.8  # phi[i][j]

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    bacterium = TB_Model.Bacterium((1, 1), 'fast')
    self.topology.automata[0].bacteria.append(bacterium)
    self.topology.automata[0].grid[1, 1]['contents'] = bacterium
    self.topology.automata[0].grid[1, 1]['oxygen'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 1]['oxygen'], 10.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[1, 1]['oxygen'], 10 - 0.001 * (4.8 * 10))


def test_chemotherapy_diffusion_no_changes(self):
    """
    No diffusion, no oxygen from source. Starts zero, stays zero
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 0.0
    self.params['chemotherapy_from_source'] = 0.0
    self.params['chemotherapy_decay'] = 0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.0)


def test_chemotherapy_not_in_chemo_window(self):
    """
    No chemo when not in window (False passed as argument to diffusion)
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 5.0
    self.params['chemotherapy_from_source'] = 10.0
    self.params['chemotherapy_decay'] = 0.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(False)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.0)


def test_chemotherapy_diffusion(self):
    """
    Only diffusion (no chemo from source). Place an amount of chemo in cell and check diffusions
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 1.0
    self.params['chemotherapy_from_source'] = 0.0
    self.params['chemotherapy_decay'] = 0.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.topology.automata[0].grid[2, 2]['chemotherapy'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 10.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 9.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.25)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.25)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.25)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.25)

    # Double diffusion rate, should double amount diffused
    self.params['chemotherapy_diffusion'] = 2.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.topology.automata[0].grid[2, 2]['chemotherapy'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 10.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 8.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.5)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.5)


def test_chemotherapy_diffusion_different_rates(self):
    """
    Only diffusion (no chemo from source). Place an amount of chemo in cell and check diffusions
    :return:
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 1.0
    self.params['chemotherapy_from_source'] = 0.0
    self.params['chemotherapy_decay'] = 0.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.topology.automata[0].grid[2, 2]['chemotherapy'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 10.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].grid[2, 1]['chemotherapy_diffusion_rate'] = 0.5
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 9.0625)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.25)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.1875)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.25)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.25)


def test_chemotherapy_from_source(self):
    """
    Only chemo is from source, no diffusion
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 0.0
    self.params['chemotherapy_from_source'] = 10.0
    self.params['chemotherapy_decay'] = 0.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 1.5 * 10.0 * 0.001)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.0)

    # Change value from source
    self.params['chemotherapy_from_source'] = 3.5

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.assertEqual(self.topology.automata[0].grid[1, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 1]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 2]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
    self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[2, 2]['chemotherapy'], 1.5 * 3.5 * 0.001)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[1, 2]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 1]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[2, 3]['chemotherapy'], 0.0)
    self.assertAlmostEqual(self.topology.automata[0].work_grid[3, 2]['chemotherapy'], 0.0)


def test_chemotherapy_decay(self):
    """
    No diffusion or from source. Initial chemo placed in cell decays
    """
    self.params['blood_vessel_value'] = 1.5  # m[i][j]
    self.params['spatial_step'] = 0.2  # dx or dy
    self.params['time_step'] = 0.001  # dt

    self.params['chemotherapy_diffusion'] = 0.0
    self.params['chemotherapy_from_source'] = 0.0
    self.params['chemotherapy_decay'] = 2.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.topology.automata[0].grid[1, 1]['chemotherapy'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 1]['chemotherapy'], 10.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[1, 1]['chemotherapy'], 10 - 0.001 * (2.0 * 10.0))

    # Change decay value
    self.params['chemotherapy_decay'] = 8.0

    self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.atts, self.params, [(2, 2)], [], [], [])
    self.sort_out_halo()

    self.topology.automata[0].grid[1, 1]['chemotherapy'] = 10.0

    self.assertEqual(self.topology.automata[0].grid[1, 1]['chemotherapy'], 10.0)

    self.topology.automata[0].diffusion_pre_process()
    self.topology.automata[0].diffusion(True)

    self.assertEqual(self.topology.automata[0].work_grid[1, 1]['chemotherapy'], 10 - 0.001 * (8.0 * 10.0))


# TODO - chemokine tests

if __name__ == '__main__':
    unittest.main()
