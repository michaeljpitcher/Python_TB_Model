import unittest
import TB_Model
import math
import numpy as np
import os


class TileTestCase(unittest.TestCase):

    def setUp(self):
        self.tile = TB_Model.Tile([5, 5], ['a', 'b', 'c'])
        self.tile.grid[0, 0]['a'] = 1.1
        self.tile.create_work_grid()

    def test_shape(self):
        self.assertEqual(self.tile.grid.shape, (5, 5))

    def test_create_work_grid(self):

        self.assertEqual(self.tile.work_grid[0, 0]['a'], 1.1)

        self.tile.grid[1, 1]['b'] = 99.9
        self.assertEqual(self.tile.work_grid[1, 1]['b'], 0.0)

    def test_attributes(self):
        self.assertTrue(isinstance(self.tile.grid[0, 0], dict))

        for x in range(5):
            for y in range(5):
                self.assertItemsEqual(self.tile.grid[x, y].keys(), ['a', 'b', 'c'])
                self.assertItemsEqual(self.tile.work_grid[x, y].keys(), ['a', 'b', 'c'])

    def test_is_on_grid(self):
        for x in range(5):
            for y in range(5):
                self.assertTrue(self.tile.address_is_on_grid([x, y]))

        for x in range(-5, 0):
            for y in range(-5, 0):
                self.assertFalse(self.tile.address_is_on_grid([x, y]))

        for x in range(5, 10):
            for y in range(5, 10):
                self.assertFalse(self.tile.address_is_on_grid([x, y]))

    def test_integer_to_address(self):

        location = 0
        for x in range(5):
            for y in range(5):
                self.assertEqual(self.tile.location_to_address(location), [x, y])
                location += 1

    def test_swap_grids(self):
        self.tile.set_attribute_work_grid([0, 1], 'a', 1)
        self.tile.swap_grids()
        self.assertEqual(self.tile.grid[0, 1]['a'], 1)
        self.assertEqual(self.tile.work_grid[0, 1]['a'], 0.0)

    def test_set_DZ_addresses(self):
        self.tile.set_addresses_for_danger_zone([[0, 0], [1, 1], [2, 2]])
        self.assertSequenceEqual(self.tile.danger_zone_addresses, [[0, 0], [1, 1], [2, 2]])

    def test_get_DZ(self):

        self.tile.set_addresses_for_danger_zone([[4, 4], [4, 3], [3, 4]])
        self.tile.grid[4, 4]['a'] = 99
        self.tile.grid[4, 3]['b'] = 99
        self.tile.grid[3, 4]['c'] = 99

        danger_zone = self.tile.get_danger_zone()

        self.assertEqual(len(danger_zone), 3)
        self.assertEqual(danger_zone[0]['a'], 99)
        self.assertEqual(danger_zone[0]['b'], 0.0)
        self.assertEqual(danger_zone[0]['c'], 0.0)
        self.assertEqual(danger_zone[1]['a'], 0.0)
        self.assertEqual(danger_zone[1]['b'], 99)
        self.assertEqual(danger_zone[1]['c'], 0.0)
        self.assertEqual(danger_zone[2]['a'], 0.0)
        self.assertEqual(danger_zone[2]['b'], 0.0)
        self.assertEqual(danger_zone[2]['c'], 99)

    def test_set_addresses_for_halo(self):
        self.tile.set_addresses_for_halo([[-1, -1], [-1, 0], [0, -1]])
        self.assertSequenceEqual(self.tile.halo_addresses, [[-1, -1], [-1, 0], [0, -1]])

    def test_set_halo(self):
        self.tile.set_addresses_for_halo([[-1, -1], [-1, 0], [0, -1]])
        cells = []
        for i in range(3):
            cell = dict()
            cell['a'] = i
            cell['b'] = i
            cell['c'] = i
            cells.append(cell)
        self.tile.set_halo(cells)
        self.assertEqual(self.tile.halo_cells[0]['a'], 0)
        self.assertEqual(self.tile.halo_cells[0]['b'], 0)
        self.assertEqual(self.tile.halo_cells[0]['c'], 0)
        self.assertEqual(self.tile.halo_cells[1]['a'], 1)
        self.assertEqual(self.tile.halo_cells[1]['b'], 1)
        self.assertEqual(self.tile.halo_cells[1]['c'], 1)
        self.assertEqual(self.tile.halo_cells[2]['a'], 2)
        self.assertEqual(self.tile.halo_cells[2]['b'], 2)
        self.assertEqual(self.tile.halo_cells[2]['c'], 2)

    def test_get_on_grid(self):
        self.assertItemsEqual(self.tile.get([0, 1]).keys(), ['a', 'b', 'c'])
        self.assertItemsEqual(self.tile.get([0, 1]).values(), [0.0, 0.0, 0.0])

    def test_get_attribute_on_grid(self):
        self.assertEqual(self.tile.get_attribute([0, 1], 'a'), 0.0)

    def test_get_on_halo(self):
        self.tile.set_addresses_for_halo([[-1, -1], [-1, 0], [0, -1]])
        cells = []
        for i in range(3):
            cell = dict()
            cell['a'] = i
            cell['b'] = i
            cell['c'] = i
            cells.append(cell)
        self.tile.set_halo(cells)

        self.assertItemsEqual(self.tile.get([-1, -1]).keys(), ['a', 'b', 'c'])

    def test_get_attribute_on_halo(self):
        self.tile.set_addresses_for_halo([[-1, -1], [-1, 0], [0, -1]])
        cells = []
        for i in range(3):
            cell = dict()
            cell['a'] = i
            cell['b'] = i
            cell['c'] = i
            cells.append(cell)
        self.tile.set_halo(cells)

        self.assertEqual(self.tile.get_attribute([-1, -1], 'a'), 0.0)

    def test_set_att(self):

        self.tile.set_attribute_grid([2, 2], 'a', 99.0)
        self.assertEqual(self.tile.grid[2, 2]['a'], 99.0)
        # Check an exception is raised if the attribute doesn't exist
        self.assertRaises(Exception, self.tile.set_attribute_grid, [2, 2], 'not set', 99.0)

    def test_set_att_work(self):
        self.tile.set_attribute_work_grid([2, 2], 'a', 99.0)
        self.assertEqual(self.tile.work_grid[2, 2]['a'], 99.0)
        # Check an exception is raised if the attribute doesn't exist
        self.assertRaises(Exception, self.tile.set_attribute_grid, [2, 2], 'not set', 99.0)


class NeighbourhoodTestCase(unittest.TestCase):

    def setUp(self):
        self.neighbourhood = TB_Model.Neighbourhood(2, 3)

    def test_neighbour_table_moore(self):
        table = self.neighbourhood.neighbour_table_moore[1]
        self.assertItemsEqual(table, [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)])

        table = self.neighbourhood.neighbour_table_moore[2]
        self.assertItemsEqual(table, [(-2, -2), (-2, -1), (-2, 0), (-2, 1), (-2, 2), (-1, -2), (-1, 2), (0, -2), (0, 2),
                                      (1, -2), (1, 2), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)])

        table = self.neighbourhood.neighbour_table_moore[3]
        self.assertItemsEqual(table,
                              [(-3, -3), (-3, -2), (-3, -1), (-3, 0), (-3, 1), (-3, 2), (-3, 3), (-2, -3), (-2, 3),
                               (-1, -3), (-1, 3), (0, -3), (0, 3), (1, -3), (1, 3), (2, -3), (2, 3), (3, -3), (3, -2),
                               (3, -1), (3, 0), (3, 1), (3, 2), (3, 3)])

    def test_calculate_neighbour_location_moore(self):
        table = self.neighbourhood.neighbour_table_moore[1]
        locations = self.neighbourhood.calculate_neighbours_locations([2, 2], table)
        self.assertItemsEqual(locations, [[1, 1], [1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2], [3, 3]])

    def test_calculate_neighbour_location_von_neumann(self):
        table = self.neighbourhood.neighbour_table_von_neumann[1]
        locations = self.neighbourhood.calculate_neighbours_locations([2, 2], table)
        self.assertItemsEqual(locations, [[2, 3], [1, 2], [2, 1], [3, 2]])

    def test_neighbours_moore(self):
        nm_2_exc = self.neighbourhood.neighbours_moore([2, 2], 2)
        self.assertItemsEqual(nm_2_exc, [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4],
                                         [1, 0], [1, 4],
                                         [2, 0], [2, 4],
                                         [3, 0], [3, 4],
                                         [4, 0], [4, 1], [4, 2], [4, 3], [4, 4]])

    def test_neighbours_von_neumann(self):
        nvn_2_exc = self.neighbourhood.neighbours_von_neumann([2, 2], 2)
        self.assertItemsEqual(nvn_2_exc, [[0, 2], [1, 1], [1, 3], [2, 0], [2, 4], [3, 1], [3, 3], [4, 2]])


class AutomatonTestCase(unittest.TestCase):

    def setUp(self):
        # TODO - more tests
        params = dict()
        params['max_depth'] = 3
        params['caseum_distance_to_reduce_diffusion'] = 1
        atts = ['a', 'contents', 'chemokine']
        self.automaton = TB_Model.Automaton([5, 5], 1, atts, params, [], [], [], [])

    def test_find_max_chemokine_neighbour(self):
        self.automaton.grid[1, 1]['chemokine'] = 99
        self.automaton.set_max_chemokine_global(99)
        neighbours = self.automaton.neighbours_moore([2, 2], 1)
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[0], neighbours.index([1, 1]))
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[1], 100)

    def test_find_max_chemokine_neighbour_tied(self):
        # [1,1] and [3,3] are equal
        self.automaton.grid[1, 1]['chemokine'] = 99
        self.automaton.grid[3, 3]['chemokine'] = 99
        self.automaton.set_max_chemokine_global(99)
        neighbours = self.automaton.neighbours_moore([2, 2], 1)
        # Force random to be [3,3]
        np.random.seed(1)
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[0], neighbours.index([3, 3]))
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[1], 100)
        # Force random to be [1,1]
        np.random.seed(100)
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[0], neighbours.index([1, 1]))
        self.assertEqual(self.automaton.find_max_chemokine_neighbour(neighbours)[1], 100)

    def test_persist_agents(self):
        b = TB_Model.Bacterium([0, 0], 'fast')
        m = TB_Model.Macrophage([0, 1], 'resting')
        t = TB_Model.TCell([0, 2])

        self.automaton.bacteria.append(b)
        self.automaton.macrophages.append(m)
        self.automaton.t_cells.append(t)

        self.automaton.persist_agents()

        self.assertEqual(self.automaton.work_grid[0, 0]['contents'], b)
        self.assertEqual(self.automaton.work_grid[0, 1]['contents'], m)
        self.assertEqual(self.automaton.work_grid[0, 2]['contents'], t)

class TopologyTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.0
        params['caseum_distance_to_reduce_diffusion'] = 0
        atts = ['a', 'b', 'c', 'blood_vessel', 'contents', 'oxygen']

        bv = [[[1, 1]], [], [], []]
        fb = [[], [[1, 1]], [], []]
        sb = [[], [], [], [[3, 3]]]
        m = [[], [], [[2, 2]], []]

        self.topology = TB_Model.Topology([2, 2], [10, 10], atts, params, bv, fb, sb, m)

    def test_init(self):
        self.assertEqual(len(self.topology.automata), 4)
        self.assertEqual(self.topology.automata[0].shape[0], 5)
        self.assertEqual(self.topology.automata[0].shape[1], 5)
        self.assertEqual(self.topology.automata[1].shape[0], 5)
        self.assertEqual(self.topology.automata[1].shape[1], 5)
        self.assertEqual(self.topology.automata[2].shape[0], 5)
        self.assertEqual(self.topology.automata[2].shape[1], 5)
        self.assertEqual(self.topology.automata[3].shape[0], 5)
        self.assertEqual(self.topology.automata[3].shape[1], 5)
        addresses = [[-3, -3], [-3, -2], [-3, -1], [-3, 0], [-3, 1], [-3, 2], [-3, 3], [-2, -3], [-2, -2], [-2, -1],
                     [-2, 0],
                     [-2, 1], [-2, 2], [-2, 3], [-1, -3], [-1, -2], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [-1, 3],
                     [0, -3], [0, -2],
                     [0, -1], [1, -3], [1, -2], [1, -1], [2, -3], [2, -2], [2, -1], [3, -3], [3, -2], [3, -1], [-3, 4],
                     [-2, 4],
                     [-1, 4], [-3, 5], [-2, 5], [-1, 5], [0, 5], [1, 5], [2, 5], [3, 5], [-3, 6], [-2, 6], [-1, 6],
                     [0, 6], [1, 6],
                     [2, 6], [3, 6], [-3, 7], [-2, 7], [-1, 7], [0, 7], [1, 7], [2, 7], [3, 7], [4, -3], [4, -2],
                     [4, -1], [4, 5],
                     [4, 6], [4, 7], [5, -3], [5, -2], [5, -1], [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5], [5, 6],
                     [5, 7],
                     [6, -3], [6, -2], [6, -1], [6, 0], [6, 1], [6, 2], [6, 3], [6, 4], [6, 5], [6, 6], [6, 7], [7, -3],
                     [7, -2],
                     [7, -1], [7, 0], [7, 1], [7, 2], [7, 3], [7, 4], [7, 5], [7, 6], [7, 7]]
        self.assertItemsEqual(self.topology.automata[0].halo_addresses, addresses)
        self.assertItemsEqual(self.topology.automata[1].halo_addresses, addresses)
        self.assertItemsEqual(self.topology.automata[2].halo_addresses, addresses)
        self.assertItemsEqual(self.topology.automata[3].halo_addresses, addresses)
        addresses = [[-1, -1], [-1, 0], [-1, 1], [-1, 2], [-1, 3],
                     [0, -1], [1, -1], [2, -1], [3, -1],
                     [-1, 4], [-1, 5], [0, 5], [1, 5], [2, 5], [3, 5],
                     [4, -1], [4, 5],
                     [5, -1], [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5]]
        self.assertItemsEqual(self.topology.automata[0].halo_depth1, addresses)
        self.assertItemsEqual(self.topology.automata[1].halo_depth1, addresses)
        self.assertItemsEqual(self.topology.automata[2].halo_depth1, addresses)
        self.assertItemsEqual(self.topology.automata[3].halo_depth1, addresses)

    def test_get_external_addresses_required(self):
        addresses = [[-1, -1], [-1, 0], [-1, 1], [-1, 2], [-1, 3],
                     [0, -1], [1, -1], [2, -1], [3, -1],
                     [-1, 4], [-1, 5], [0, 5], [1, 5], [2, 5], [3, 5],
                     [4, -1], [4, 5],
                     [5, -1], [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5]
                     ]
        gear = self.topology.get_external_addresses_required(self.topology.automata[0], 1)
        self.assertItemsEqual(addresses, gear)

        addresses = [[-3, -3], [-3, -2], [-3, -1], [-3, 0], [-3, 1], [-3, 2], [-3, 3], [-2, -3], [-2, -2], [-2, -1],
                     [-2, 0],
                     [-2, 1], [-2, 2], [-2, 3], [-1, -3], [-1, -2], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [-1, 3],
                     [0, -3], [0, -2],
                     [0, -1], [1, -3], [1, -2], [1, -1], [2, -3], [2, -2], [2, -1], [3, -3], [3, -2], [3, -1], [-3, 4],
                     [-2, 4],
                     [-1, 4], [-3, 5], [-2, 5], [-1, 5], [0, 5], [1, 5], [2, 5], [3, 5], [-3, 6], [-2, 6], [-1, 6],
                     [0, 6], [1, 6],
                     [2, 6], [3, 6], [-3, 7], [-2, 7], [-1, 7], [0, 7], [1, 7], [2, 7], [3, 7], [4, -3], [4, -2],
                     [4, -1], [4, 5],
                     [4, 6], [4, 7], [5, -3], [5, -2], [5, -1], [5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5], [5, 6],
                     [5, 7],
                     [6, -3], [6, -2], [6, -1], [6, 0], [6, 1], [6, 2], [6, 3], [6, 4], [6, 5], [6, 6], [6, 7], [7, -3],
                     [7, -2],
                     [7, -1], [7, 0], [7, 1], [7, 2], [7, 3], [7, 4], [7, 5], [7, 6], [7, 7]]

        gear = self.topology.get_external_addresses_required(self.topology.automata[0], 3)
        self.assertItemsEqual(addresses, gear)


class TwoDimensionalTopologyTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.0
        params['caseum_distance_to_reduce_diffusion'] = 0
        atts = ['a', 'b', 'c', 'blood_vessel', 'contents', 'oxygen']
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, [[3, 3]], [[1, 1]], [[9, 9]],
                                                        [[7, 1]])

    def test_init(self):
        self.assertSequenceEqual(self.topology.origins, [[0, 0], [0, 5], [5, 0], [5, 5]])

        tile0_dz_addresses = [[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4], [2, 0], [2, 1], [2, 2], [2, 3], [2, 4],
                              [3, 0], [3, 1], [3, 2], [3, 3], [3, 4], [4, 0], [4, 1], [4, 2], [4, 3], [4, 4]]

        self.assertItemsEqual(tile0_dz_addresses, self.topology.automata[0].danger_zone_addresses)
        tile1_dz_addresses = [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2], [2, 3], [2, 4],
                              [3, 0], [3, 1], [3, 2], [3, 3], [3, 4], [4, 0], [4, 1], [4, 2], [4, 3], [4, 4]]
        self.assertItemsEqual(tile1_dz_addresses, self.topology.automata[1].danger_zone_addresses)
        tile2_dz_addresses = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [2, 0],
                              [2, 1], [2, 2], [2, 3], [2, 4], [3, 2], [3, 3], [3, 4], [4, 2], [4, 3], [4, 4]]
        self.assertItemsEqual(tile2_dz_addresses, self.topology.automata[2].danger_zone_addresses)
        tile3_dz_addresses = [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [2, 0],
                              [2, 1], [2, 2], [2, 3], [2, 4], [3, 0], [3, 1], [3, 2], [4, 0], [4, 1], [4, 2]]
        self.assertItemsEqual(tile3_dz_addresses, self.topology.automata[3].danger_zone_addresses)

    def test_normalise_address(self):
        self.assertSequenceEqual(self.topology.normalise_address([0, 0]), [0, 0])
        self.assertSequenceEqual(self.topology.normalise_address([9, 9]), [9, 9])
        self.assertEqual(self.topology.normalise_address([-1, 0]), None)
        self.assertEqual(self.topology.normalise_address([3, 10]), None)

    def test_blood_vessel(self):

        for a in self.topology.automata:
            for x in range(5):
                for y in range(5):
                    if a.tile_id == 0 and x == 3 and y == 3:
                        self.assertEqual(a.grid[x, y]['blood_vessel'], 1.5)
                        self.assertEqual(a.work_grid[x, y]['blood_vessel'], 1.5)
                    else:
                        self.assertEqual(a.grid[x, y]['blood_vessel'], 0.0)
                        self.assertEqual(a.work_grid[x, y]['blood_vessel'], 0.0)

    def test_global_to_local(self):
        self.assertEqual(self.topology.global_to_local([0, 0])[0], 0)
        self.assertSequenceEqual(self.topology.global_to_local([0, 0])[1], [0, 0])
        self.assertEqual(self.topology.global_to_local([2, 6])[0], 1)
        self.assertSequenceEqual(self.topology.global_to_local([2, 6])[1], [2, 1])
        self.assertEqual(self.topology.global_to_local([7, 3])[0], 2)
        self.assertSequenceEqual(self.topology.global_to_local([7, 3])[1], [2, 3])
        self.assertEqual(self.topology.global_to_local([8, 8])[0], 3)
        self.assertSequenceEqual(self.topology.global_to_local([8, 8])[1], [3, 3])

    def test_local_to_global(self):
        self.assertEqual(self.topology.local_to_global(0, [-1, -1]), None)
        self.assertEqual(self.topology.local_to_global(1, [-1, 0]), None)
        self.assertEqual(self.topology.local_to_global(1, [2, 7]), None)
        self.assertEqual(self.topology.local_to_global(2, [0, -1]), None)
        self.assertEqual(self.topology.local_to_global(2, [8, 3]), None)
        self.assertEqual(self.topology.local_to_global(3, [0, 8]), None)
        self.assertEqual(self.topology.local_to_global(3, [8, 1]), None)
        self.assertSequenceEqual(self.topology.local_to_global(1, [0, -1]), [0, 4])
        self.assertSequenceEqual(self.topology.local_to_global(1, [5, 0]), [5, 5])
        self.assertSequenceEqual(self.topology.local_to_global(1, [7, -2]), [7, 3])
        self.assertSequenceEqual(self.topology.local_to_global(2, [-2, 2]), [3, 2])
        self.assertSequenceEqual(self.topology.local_to_global(2, [0, 6]), [5, 6])
        self.assertSequenceEqual(self.topology.local_to_global(2, [-1, 6]), [4, 6])
        self.assertSequenceEqual(self.topology.local_to_global(3, [-1, -1]), [4, 4])
        self.assertSequenceEqual(self.topology.local_to_global(3, [-2, 3]), [3, 8])
        self.assertSequenceEqual(self.topology.local_to_global(3, [2, -3]), [7, 2])

    def test_local_to_local(self):
        self.attributes = ['a', 'blood_vessel', 'oxygen']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.parameters['initial_oxygen'] = 1.0
        self.parameters['caseum_distance_to_reduce_diffusion'] = 0
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters, [[3, 3]], [],
                                                        [], [])

        self.assertSequenceEqual(self.topology.local_to_local(0, [0, 4], 1), [0, 1])
        self.assertSequenceEqual(self.topology.local_to_local(0, [3, 2], 2), [0, 2])
        self.assertSequenceEqual(self.topology.local_to_local(0, [3, 3], 3), [0, 0])

        self.assertSequenceEqual(self.topology.local_to_local(0, [0, 2], 1), [0, -1])
        self.assertSequenceEqual(self.topology.local_to_local(1, [0, 0], 2), [-3, 3])

        self.assertSequenceEqual(self.topology.local_to_local(3, [-1, -1], 0), [2, 2])
        self.assertSequenceEqual(self.topology.local_to_local(3, [-1, 0], 1), [2, 0])
        self.assertSequenceEqual(self.topology.local_to_local(3, [0, -2], 2), [0, 1])

    def test_create_halos(self):
        self.attributes = ['a', 'blood_vessel', 'oxygen']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.parameters['initial_oxygen'] = 1.0
        self.parameters['caseum_distance_to_reduce_diffusion'] = 0
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters, [[3, 3]], [],
                                                        [], [])

        tile_id = 0
        x = 0
        y = 0
        for i in range(36):
            self.topology.automata[tile_id].grid[x, y]['a'] = i
            y += 1
            if y == 3 and x == 2:
                tile_id += 1
                y = x = 0
            elif y == 3:
                y = 0
                x += 1

        danger_zone_values = []
        for a in range(4):
            danger_zone_values.append(self.topology.automata[a].get_danger_zone())

        halos = self.topology.create_halos(danger_zone_values)

        for index in range(4):
            self.assertEqual(len(halos[index]), 16)
            if index == 0:
                addresses = self.topology.automata[index].halo_addresses
                self.assertEqual(halos[index][addresses.index([0, 3])]['a'], 9)
                self.assertEqual(halos[index][addresses.index([1, 3])]['a'], 12)
                self.assertEqual(halos[index][addresses.index([2, 3])]['a'], 15)
                self.assertEqual(halos[index][addresses.index([3, 0])]['a'], 18)
                self.assertEqual(halos[index][addresses.index([3, 1])]['a'], 19)
                self.assertEqual(halos[index][addresses.index([3, 2])]['a'], 20)
                self.assertEqual(halos[index][addresses.index([3, 3])]['a'], 27)
                self.assertEqual(halos[index][addresses.index([-1, -1])], None)
                self.assertEqual(halos[index][addresses.index([-1, 0])], None)
                self.assertEqual(halos[index][addresses.index([-1, 1])], None)
                self.assertEqual(halos[index][addresses.index([-1, 2])], None)
                self.assertEqual(halos[index][addresses.index([-1, 3])], None)
                self.assertEqual(halos[index][addresses.index([0, -1])], None)
                self.assertEqual(halos[index][addresses.index([1, -1])], None)
                self.assertEqual(halos[index][addresses.index([2, -1])], None)
                self.assertEqual(halos[index][addresses.index([3, -1])], None)
            elif index == 1:
                self.assertEqual(halos[index][addresses.index([0, 3])], None)
                self.assertEqual(halos[index][addresses.index([1, 3])], None)
                self.assertEqual(halos[index][addresses.index([2, 3])], None)
                self.assertEqual(halos[index][addresses.index([3, 0])]['a'], 27)
                self.assertEqual(halos[index][addresses.index([3, 1])]['a'], 28)
                self.assertEqual(halos[index][addresses.index([3, 2])]['a'], 29)
                self.assertEqual(halos[index][addresses.index([3, 3])], None)
                self.assertEqual(halos[index][addresses.index([-1, -1])], None)
                self.assertEqual(halos[index][addresses.index([-1, 0])], None)
                self.assertEqual(halos[index][addresses.index([-1, 1])], None)
                self.assertEqual(halos[index][addresses.index([-1, 2])], None)
                self.assertEqual(halos[index][addresses.index([-1, 3])], None)
                self.assertEqual(halos[index][addresses.index([0, -1])]['a'], 2)
                self.assertEqual(halos[index][addresses.index([1, -1])]['a'], 5)
                self.assertEqual(halos[index][addresses.index([2, -1])]['a'], 8)
                self.assertEqual(halos[index][addresses.index([3, -1])]['a'], 20)
            elif index == 2:
                self.assertEqual(halos[index][addresses.index([0, 3])]['a'], 27)
                self.assertEqual(halos[index][addresses.index([1, 3])]['a'], 30)
                self.assertEqual(halos[index][addresses.index([2, 3])]['a'], 33)
                self.assertEqual(halos[index][addresses.index([3, 0])], None)
                self.assertEqual(halos[index][addresses.index([3, 1])], None)
                self.assertEqual(halos[index][addresses.index([3, 2])], None)
                self.assertEqual(halos[index][addresses.index([3, 3])], None)
                self.assertEqual(halos[index][addresses.index([-1, -1])], None)
                self.assertEqual(halos[index][addresses.index([-1, 0])]['a'], 6)
                self.assertEqual(halos[index][addresses.index([-1, 1])]['a'], 7)
                self.assertEqual(halos[index][addresses.index([-1, 2])]['a'], 8)
                self.assertEqual(halos[index][addresses.index([-1, 3])]['a'], 15)
                self.assertEqual(halos[index][addresses.index([0, -1])], None)
                self.assertEqual(halos[index][addresses.index([1, -1])], None)
                self.assertEqual(halos[index][addresses.index([2, -1])], None)
                self.assertEqual(halos[index][addresses.index([3, -1])], None)
            elif index == 3:
                self.assertEqual(halos[index][addresses.index([0, 3])], None)
                self.assertEqual(halos[index][addresses.index([1, 3])], None)
                self.assertEqual(halos[index][addresses.index([2, 3])], None)
                self.assertEqual(halos[index][addresses.index([3, 0])], None)
                self.assertEqual(halos[index][addresses.index([3, 1])], None)
                self.assertEqual(halos[index][addresses.index([3, 2])], None)
                self.assertEqual(halos[index][addresses.index([3, 3])], None)
                self.assertEqual(halos[index][addresses.index([-1, -1])]['a'], 8)
                self.assertEqual(halos[index][addresses.index([-1, 0])]['a'], 15)
                self.assertEqual(halos[index][addresses.index([-1, 1])]['a'], 16)
                self.assertEqual(halos[index][addresses.index([-1, 2])]['a'], 17)
                self.assertEqual(halos[index][addresses.index([-1, 3])], None)
                self.assertEqual(halos[index][addresses.index([0, -1])]['a'], 20)
                self.assertEqual(halos[index][addresses.index([1, -1])]['a'], 23)
                self.assertEqual(halos[index][addresses.index([2, -1])]['a'], 26)
                self.assertEqual(halos[index][addresses.index([3, -1])], None)

    def test_get_external(self):
        self.attributes = ['a', 'blood_vessel', 'oxygen']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.parameters['initial_oxygen'] = 1.0
        self.parameters['caseum_distance_to_reduce_diffusion'] = 0
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters, [[3, 3]], [],
                                                        [], [])

        tile_id = 0
        x = 0
        y = 0
        for i in range(36):
            self.topology.automata[tile_id].grid[x, y]['a'] = i
            y += 1
            if y == 3 and x == 2:
                tile_id += 1
                y = x = 0
            elif y == 3:
                y = 0
                x += 1

        danger_zone_values = []
        for a in range(4):
            danger_zone_values.append(self.topology.automata[a].get_danger_zone())

        halos = self.topology.create_halos(danger_zone_values)

        for index in range(4):
            self.topology.automata[index].set_halo(halos[index])
            if index == 0:
                self.assertEqual(self.topology.automata[index].get_attribute([0, 3], 'a'), 9)
                self.assertEqual(self.topology.automata[index].get_attribute([1, 3], 'a'), 12)
                self.assertEqual(self.topology.automata[index].get_attribute([2, 3], 'a'), 15)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 0], 'a'), 18)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 1], 'a'), 19)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 2], 'a'), 20)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 3], 'a'), 27)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 0], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 2], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([0, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([1, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([2, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, -1], 'a'), None)
            elif index == 1:
                self.assertEqual(self.topology.automata[index].get_attribute([0, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([1, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([2, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 0], 'a'), 27)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 1], 'a'), 28)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 2], 'a'), 29)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 0], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 2], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([0, -1], 'a'), 2)
                self.assertEqual(self.topology.automata[index].get_attribute([1, -1], 'a'), 5)
                self.assertEqual(self.topology.automata[index].get_attribute([2, -1], 'a'), 8)
                self.assertEqual(self.topology.automata[index].get_attribute([3, -1], 'a'), 20)

            elif index == 2:
                self.assertEqual(self.topology.automata[index].get_attribute([0, 3], 'a'), 27)
                self.assertEqual(self.topology.automata[index].get_attribute([1, 3], 'a'), 30)
                self.assertEqual(self.topology.automata[index].get_attribute([2, 3], 'a'), 33)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 0], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 2], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 0], 'a'), 6)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 1], 'a'), 7)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 2], 'a'), 8)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 3], 'a'), 15)
                self.assertEqual(self.topology.automata[index].get_attribute([0, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([1, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([2, -1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, -1], 'a'), None)

            elif index == 3:
                self.assertEqual(self.topology.automata[index].get_attribute([0, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([1, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([2, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 0], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 1], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 2], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([3, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, -1], 'a'), 8)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 0], 'a'), 15)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 1], 'a'), 16)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 2], 'a'), 17)
                self.assertEqual(self.topology.automata[index].get_attribute([-1, 3], 'a'), None)
                self.assertEqual(self.topology.automata[index].get_attribute([0, -1], 'a'), 20)
                self.assertEqual(self.topology.automata[index].get_attribute([1, -1], 'a'), 23)
                self.assertEqual(self.topology.automata[index].get_attribute([2, -1], 'a'), 26)
                self.assertEqual(self.topology.automata[index].get_attribute([3, -1], 'a'), None)


class TBAutomatonScenariosTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['caseum_threshold_to_reduce_diffusion'] = 2
        params['oxygen_diffusion_caseum_reduction'] = 1.5
        params['chemotherapy_diffusion_caseum_reduction'] = 1.5
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['oxygen_uptake_from_bacteria'] = 1.0
        params['time_step'] = 0.001
        params['chemotherapy_from_source'] = 1.0
        params['chemotherapy_decay'] = 0.35
        params['chemokine_diffusion'] = 0.05
        params['chemokine_from_bacteria'] = 0.5
        params['chemokine_from_macrophages'] = 1
        params['chemokine_decay'] = 0.347

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[3, 3]]
        fast_bacteria = [[1, 1]]
        slow_bacteria = [[9, 9]]
        macrophages = [[7, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

        # Create a halo - not needed for tests but needed for code
        halo = []
        for a in self.topology.automata[0].halo_addresses:
            x, y = a
            if x < 0 or y < 0:
                halo.append(None)
            else:
                cell = dict()
                cell['blood_vessel'] = 0.0
                cell['contents'] = 0.0
                cell['oxygen'] = 0.0
                cell['chemotherapy'] = 0.0
                cell['chemokine'] = 0.0
                cell['oxygen_diffusion_rate'] = 0.0
                cell['chemotherapy_diffusion_rate'] = 0.0
                halo.append(cell)
        self.topology.automata[0].set_halo(halo)

    def test_initialise(self):
        self.assertEqual(len(self.topology.automata[0].blood_vessels), 1)
        self.assertItemsEqual(self.topology.automata[0].blood_vessels[0], [3, 3])
        self.assertEqual(len(self.topology.automata[1].blood_vessels), 0)
        self.assertEqual(len(self.topology.automata[2].blood_vessels), 0)
        self.assertEqual(len(self.topology.automata[3].blood_vessels), 0)

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)
        self.assertEqual(self.topology.automata[0].bacteria[0].metabolism, "fast")
        self.assertItemsEqual(self.topology.automata[0].bacteria[0].address, [1, 1])
        self.assertEqual(len(self.topology.automata[1].bacteria), 0)
        self.assertEqual(len(self.topology.automata[2].bacteria), 0)
        self.assertEqual(len(self.topology.automata[3].bacteria), 1)
        self.assertEqual(self.topology.automata[3].bacteria[0].metabolism, "slow")
        self.assertItemsEqual(self.topology.automata[3].bacteria[0].address, [4, 4])

    def test_pre_process_caseum(self):

        # Add some caseum to automaton 0
        self.topology.automata[0].grid[0, 0]['contents'] = 'caseum'
        self.topology.automata[0].caseum.append([0, 0])
        self.topology.automata[0].grid[0, 1]['contents'] = 'caseum'
        self.topology.automata[0].caseum.append([0, 1])

        # Run the pre process loop
        self.topology.automata[0].diffusion_pre_process()

        # Cells close to caseum
        self.assertEqual(self.topology.automata[0].grid[0, 2]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[1, 0]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[1, 2]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[2, 0]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[2, 1]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(self.topology.automata[0].grid[2, 2]['oxygen_diffusion_rate'], 1.0 / 1.5)

        # Cells far enough away from caseum
        self.assertEqual(self.topology.automata[0].grid[0, 3]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[0, 4]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[1, 3]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[1, 4]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[2, 4]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[3, 0]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[3, 1]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[3, 3]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(self.topology.automata[0].grid[3, 4]['oxygen_diffusion_rate'], 1.0)

    def test_pre_process_caseum_halo(self):

        # Add some caseum to automaton 0
        self.topology.automata[0].grid[0, 0]['contents'] = 'caseum'
        self.topology.automata[0].caseum.append([0, 0])
        self.topology.automata[0].grid[0, 1]['contents'] = 'caseum'
        self.topology.automata[0].caseum.append([0, 1])

        # Create a halo with caseum in [0,7] and [1,7] - which is [0,2] and [1,2] of right neighbour tile
        halo = []
        for a in self.topology.automata[0].halo_addresses:
            x, y = a
            if x < 0 or y < 0:
                halo.append(None)
            else:
                cell = dict()
                cell['blood_vessel'] = 0.0
                if (x == 0 and y == 7) or (x == 1 and y == 7):
                    cell['contents'] = 'caseum'
                else:
                    cell['contents'] = 0.0
                cell['oxygen'] = 0.0
                halo.append(cell)
        self.topology.automata[0].set_halo(halo)

        # Run the pre process loop
        self.topology.automata[0].diffusion_pre_process()

        halo_addresses = self.topology.automata[0].halo_addresses
        halo_cells = self.topology.automata[0].halo_cells

        self.assertEqual(halo_cells[halo_addresses.index([0, 5])]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(halo_cells[halo_addresses.index([1, 5])]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(halo_cells[halo_addresses.index([2, 5])]['oxygen_diffusion_rate'], 1.0 / 1.5)
        self.assertEqual(halo_cells[halo_addresses.index([3, 5])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([4, 5])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 0])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 1])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 2])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 3])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 4])]['oxygen_diffusion_rate'], 1.0)
        self.assertEqual(halo_cells[halo_addresses.index([5, 5])]['oxygen_diffusion_rate'], 1.0)

    def test_oxygen_basic(self):

        self.assertEqual(self.topology.automata[0].grid[3, 3]['oxygen'], 1.5)
        self.assertEqual(self.topology.automata[0].grid[2, 3]['oxygen'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[4, 3]['oxygen'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[3, 2]['oxygen'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[3, 4]['oxygen'], 0.0)

        self.topology.automata[0].diffusion_pre_process()

        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].oxygen(self.topology.automata[0].get([3, 3], 'grid'), neighbours), 1.3536)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([2, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].oxygen(self.topology.automata[0].get([2, 3], 'grid'), neighbours), 0.0375)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([4, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].oxygen(self.topology.automata[0].get([4, 3], 'grid'), neighbours), 0.0375)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 2], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].oxygen(self.topology.automata[0].get([3, 2], 'grid'), neighbours), 0.0375)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 4], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].oxygen(self.topology.automata[0].get([3, 4], 'grid'), neighbours), 0.0375)

    def test_chemotherapy_basic(self):

        self.assertEqual(self.topology.automata[0].grid[3, 3]['chemotherapy'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[2, 3]['chemotherapy'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[4, 3]['chemotherapy'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[3, 2]['chemotherapy'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[3, 4]['chemotherapy'], 0.0)
        self.topology.automata[0].diffusion_pre_process()

        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3,3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 3],'grid'),neighbours), 0.0015)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([2, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([2, 3],'grid'), neighbours), 0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([4, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([4, 3],'grid'),neighbours), 0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 2], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 2],'grid'),neighbours), 0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 4], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 4],'grid'),neighbours), 0.0)

        # Set value direct to grid to save time
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.topology.automata[0].set_attribute_grid([3, 3], 'chemotherapy',
                                                     self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 3], 'grid'), neighbours))

        self.topology.automata[0].diffusion_pre_process()

        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 3],'grid'),neighbours), 0.002886975)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([2, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([2, 3],'grid'),neighbours), 0.000028125)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([4, 3], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([4, 3],'grid'),neighbours), 0.000028125)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 2], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 2],'grid'),neighbours), 0.000028125)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([3, 4], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemotherapy(self.topology.automata[0].get([3, 4],'grid'),neighbours), 0.000028125)

    def test_chemokine_basic(self):
        self.topology.automata[0].diffusion_pre_process()
        self.assertEqual(self.topology.automata[0].grid[1, 1]['chemokine'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[0, 1]['chemokine'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[2, 1]['chemokine'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[1, 0]['chemokine'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[1, 2]['chemokine'], 0.0)

        self.topology.automata[0].diffusion_pre_process()
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1,1],1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 1], 'grid'), neighbours),
                         0.0005)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([0, 1], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([0, 1], 'grid'), neighbours),
                         0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([2, 1], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([2, 1], 'grid'), neighbours),
                         0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1, 0], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 0], 'grid'), neighbours),
                         0.0)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1, 2], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 2], 'grid'), neighbours),
                         0.0)

        self.topology.automata[0].set_attribute_grid([1, 1], 'chemokine', 0.0005)

        self.topology.automata[0].diffusion_pre_process()
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1, 1], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertEqual(self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 1], 'grid'), neighbours),
                         0.0009973265)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([0, 1], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].chemokine(self.topology.automata[0].get([0, 1], 'grid'), neighbours), 0.000000625)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([2, 1], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].chemokine(self.topology.automata[0].get([2, 1], 'grid'), neighbours), 0.000000625)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1, 0], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 0], 'grid'), neighbours), 0.000000625)
        neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([1, 2], 1)
        neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
        self.assertAlmostEqual(
            self.topology.automata[0].chemokine(self.topology.automata[0].get([1, 2], 'grid'), neighbours), 0.000000625)

    def test_local_and_global_levels(self):

        self.assertAlmostEqual(self.topology.automata[0].max_oxygen_local, 1.5)

        self.topology.automata[0].diffusion_pre_process()

        self.topology.automata[0].max_oxygen_local = 0.0
        self.topology.automata[0].max_chemotherapy_local = 0.0
        self.topology.automata[0].max_chemokine_local = 0.0

        for x in range(5):
            for y in range(5):
                neighbour_addresses = self.topology.automata[0].neighbours_von_neumann([x, y], 1)
                neighbours = [self.topology.automata[0].get(n) for n in neighbour_addresses]
                self.topology.automata[0].set_attribute_work_grid([x, y], 'oxygen', self.topology.automata[0].oxygen(
                    self.topology.automata[0].get([x, y], 'grid'), neighbours))
                self.topology.automata[0].set_attribute_work_grid([x, y], 'chemotherapy',
                                                                  self.topology.automata[0].chemotherapy(
                                                                      self.topology.automata[0].get([x, y], 'grid'),
                                                                      neighbours))
                self.topology.automata[0].set_attribute_work_grid([x, y], 'chemokine',
                                                                  self.topology.automata[0].chemokine(
                                                                      self.topology.automata[0].get([x, y], 'grid'),
                                                                      neighbours))

        self.topology.automata[0].swap_grids()

        self.assertAlmostEqual(self.topology.automata[0].max_oxygen_local, 1.3536)
        self.assertAlmostEqual(self.topology.automata[0].max_chemotherapy_local, 0.0015)
        self.assertAlmostEqual(self.topology.automata[0].max_chemokine_local, 0.0005)

        for a in self.topology.automata:
            a.set_max_oxygen_global(1.59)
            a.set_max_chemotherapy_global(0.5)
            a.set_max_chemokine_global(0.12)

        self.assertEqual(self.topology.automata[0].max_oxygen_global, 1.59)
        self.assertEqual(self.topology.automata[1].max_oxygen_global, 1.59)
        self.assertEqual(self.topology.automata[2].max_oxygen_global, 1.59)
        self.assertEqual(self.topology.automata[3].max_oxygen_global, 1.59)
        self.assertEqual(self.topology.automata[0].max_chemotherapy_global, 0.5)
        self.assertEqual(self.topology.automata[1].max_chemotherapy_global, 0.5)
        self.assertEqual(self.topology.automata[2].max_chemotherapy_global, 0.5)
        self.assertEqual(self.topology.automata[3].max_chemotherapy_global, 0.5)
        self.assertEqual(self.topology.automata[0].max_chemokine_global, 0.12)
        self.assertEqual(self.topology.automata[1].max_chemokine_global, 0.12)
        self.assertEqual(self.topology.automata[2].max_chemokine_global, 0.12)
        self.assertEqual(self.topology.automata[3].max_chemokine_global, 0.12)

        self.assertAlmostEqual(self.topology.automata[0].oxygen_scale([3, 3]), 1.3536/1.59 * 100)
        self.assertAlmostEqual(self.topology.automata[0].chemotherapy_scale([3, 3]), 0.0015 / 0.5 * 100)
        self.assertAlmostEqual(self.topology.automata[0].chemokine_scale([1, 1]), 0.0005 / 0.12 * 100)


# EVENT TESTING

class BacteriaReplicationTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['caseum_threshold_to_reduce_diffusion'] = 2
        params['oxygen_diffusion_caseum_reduction'] = 1.5
        params['chemotherapy_diffusion_caseum_reduction'] = 1.5
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['oxygen_uptake_from_bacteria'] = 1.0
        params['time_step'] = 1
        params['chemotherapy_from_source'] = 1.0
        params['chemotherapy_decay'] = 0.35
        params['chemokine_diffusion'] = 0.05
        params['chemokine_from_bacteria'] = 0.5
        params['chemokine_from_macrophages'] = 1
        params['chemokine_decay'] = 0.347
        params['chemotherapy_schedule1_start'] = 100
        params['chemotherapy_schedule2_start'] = 200
        params['interval_to_record_results'] = 1000

        # Set the limits to 1 and 2 - forces random number to be 1 and thus always replicates each step
        params['bacteria_replication_fast_upper'] = 2
        params['bacteria_replication_fast_lower'] = 1
        params['bacteria_replication_slow_upper'] = 2
        params['bacteria_replication_slow_lower'] = 1

        params['bacteria_threshold_for_t_cells'] = 100
        params['macrophage_recruitment_probability'] = 0
        params['chemotherapy_scale_for_kill_fast_bacteria'] = 100
        params['chemotherapy_scale_for_kill_slow_bacteria'] = 100

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[3, 3]]
        fast_bacteria = [[0, 0], [3, 5]]
        slow_bacteria = [[9, 9]]
        macrophages = [[7, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_bacteria_replication_event_creation(self):
        self.sort_out_halos()

        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumReplication))
        self.assertEqual(len(event.dependant_addresses), 2)
        self.assertSequenceEqual(event.dependant_addresses[0], [0, 0])
        self.assertTrue(event.dependant_addresses[1] == [0, 1] or event.dependant_addresses[1] == [1, 0] or
                        event.dependant_addresses[1] == [1, 1])

    def test_bacteria_replication_event_creation_slow(self):
        self.sort_out_halos()

        self.topology.automata[3].update()

        self.assertEqual(len(self.topology.automata[3].potential_events), 1)
        event = self.topology.automata[3].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumReplication))
        self.assertEqual(len(event.dependant_addresses), 2)
        self.assertSequenceEqual(event.dependant_addresses[0], [4, 4])
        self.assertTrue(event.dependant_addresses[1] == [3, 4] or event.dependant_addresses[1] == [4, 3] or
                        event.dependant_addresses[1] == [3, 3])

    def test_bacteria_replication_event_process(self):
        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumReplication))

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 2)
        self.assertTrue(isinstance(self.topology.automata[0].get_attribute(event.dependant_addresses[0], 'contents'),
                                   TB_Model.Bacterium))

        self.assertEqual(self.topology.automata[0].grid[0, 0]['contents'].division_neighbourhood, 'vn')
        self.assertEqual(self.topology.automata[0].grid[0, 0]['contents'].age, 0.0)

    def test_bacteria_replication_across_boundary(self):

        self.topology.automata[0].bacteria = []
        self.topology.automata[0].grid[0, 0]['contents'] = 0.0
        self.topology.automata[0].work_grid[0, 0]['contents'] = 0.0

        #    C C C
        #      B C
        #    C C C

        # Forces bacteria to replicate into (3,4) (globally) the only free cell in its immediate neighbourhood
        self.topology.automata[0].grid[2, 4]['contents'] = 'caseum'
        self.topology.automata[0].grid[4, 4]['contents'] = 'caseum'
        self.topology.automata[1].grid[2, 0]['contents'] = 'caseum'
        self.topology.automata[1].grid[2, 1]['contents'] = 'caseum'
        self.topology.automata[1].grid[3, 1]['contents'] = 'caseum'
        self.topology.automata[1].grid[4, 0]['contents'] = 'caseum'
        self.topology.automata[1].grid[4, 1]['contents'] = 'caseum'
        self.topology.automata[0].caseum += [[2, 4], [4, 4]]
        self.topology.automata[1].caseum += [[2, 0], [2, 1], [3, 1], [4, 0], [4, 1]]

        self.sort_out_halos()

        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(len(self.topology.automata[1].bacteria), 1)

        self.topology.automata[0].update()
        self.topology.automata[1].update()

        self.assertEqual(len(self.topology.automata[1].potential_events), 1)
        event = self.topology.automata[1].potential_events[0]
        self.assertTrue(event.dependant_addresses[0] == [3, 0])
        self.assertTrue(event.dependant_addresses[1] == [3, -1])

        new_addresses = [self.topology.local_to_local(1, event.dependant_addresses[0], 0), self.topology.local_to_local(1, event.dependant_addresses[1],0)]
        self.assertTrue(new_addresses[0] == [3, 5])
        self.assertTrue(new_addresses[1] == [3, 4])

        new_event = event.clone(new_addresses)

        self.topology.automata[0].process_events([new_event])
        self.topology.automata[1].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)
        self.assertEqual(len(self.topology.automata[1].bacteria), 1)

        self.assertTrue(isinstance(self.topology.automata[0].grid[3, 4]['contents'], TB_Model.Bacterium))
        self.assertTrue(isinstance(self.topology.automata[1].grid[3, 0]['contents'], TB_Model.Bacterium))

    def test_bacteria_replication_negative_resting(self):
        self.sort_out_halos()

        self.topology.automata[0].bacteria[0].resting = True

        for x in range(0,5):
            for y in range(0,5):
                if x != 0 or y != 0:
                    self.topology.automata[0].grid[x,y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x,y])

        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_bacteria_replication_negative_fast_age(self):
        self.sort_out_halos()

        self.topology.automata[0].parameters['bacteria_replication_fast_upper'] = 999
        self.topology.automata[0].parameters['bacteria_replication_fast_lower'] = 998

        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_bacteria_replication_negative_slow_age(self):
        self.sort_out_halos()

        self.topology.automata[3].parameters['bacteria_replication_slow_upper'] = 999
        self.topology.automata[3].parameters['bacteria_replication_slow_lower'] = 998

        self.topology.automata[3].update()

        self.assertEqual(len(self.topology.automata[3].potential_events), 0)



class TCellRecruitmentTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['time_step'] = 0.001
        params['chemokine_diffusion'] = 0.05
        params['chemokine_decay'] = 0.347
        params['chemotherapy_schedule1_start'] = 100
        params['chemotherapy_schedule2_start'] = 200
        params['macrophage_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['interval_to_record_results'] = 1000

        # PARAMETERS WHICH ARE RELEVANT FOR THESE TESTS
        # Each of these require negative testing as well (later)

        # Threshold = 0: not dependent on bacteria
        params['bacteria_threshold_for_t_cells'] = 0
        # Probability = 100: will always recruit based on random number
        params['t_cell_recruitment_probability'] = 100
        # Scale = -1, will always place in a cell regardless of scale
        params['chemokine_scale_for_t_cell_recruitment'] = -1

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[3, 3]]
        fast_bacteria = [[9, 9], [8, 9]]
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        bacteria_global = 0
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
            bacteria_global += len(i.bacteria)
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])
            self.topology.automata[i].set_global_bacteria_number(bacteria_global)

    def test_t_cell_recruited_internally(self):
        self.sort_out_halos()
        for i in range(4):
            self.assertEqual(self.topology.automata[1].number_of_bacteria_global, 2)

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        self.assertTrue(isinstance(self.topology.automata[0].potential_events[0], TB_Model.RecruitTCell))

        address = self.topology.automata[0].potential_events[0].dependant_addresses[0]
        distance = math.fabs(address[0] - 3) + math.fabs(address[1] - 3)
        self.assertEqual(distance, 1)

    def test_t_cell_recruit_process_event(self):
        self.sort_out_halos()
        self.topology.automata[0].update()
        event = self.topology.automata[0].potential_events[0]
        self.topology.automata[0].process_events([event])

        for x in range(5):
            for y in range(5):
                if x == event.dependant_addresses[0][0] and y == event.dependant_addresses[0][1]:
                    self.assertTrue(
                        isinstance(self.topology.automata[0].get_attribute([x, y], 'contents'), TB_Model.TCell))
                else:
                    self.assertEqual(self.topology.automata[0].get_attribute([x, y], 'contents'), 0.0)

    def test_t_cell_recruited_across_boundary(self):

        # Need a different set-up
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['caseum_threshold_to_reduce_diffusion'] = 200
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['time_step'] = 0.001
        params['chemokine_diffusion'] = 0.05
        params['chemokine_decay'] = 0.347
        params['chemotherapy_schedule1_start'] = 100
        params['chemotherapy_schedule2_start'] = 200
        params['macrophage_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['interval_to_record_results'] = 1000

        # PARAMETERS WHICH ARE RELEVANT FOR THESE TESTS
        params['bacteria_threshold_for_t_cells'] = 0
        params['t_cell_recruitment_probability'] = 100
        params['chemokine_scale_for_t_cell_recruitment'] = -1

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        # Only 1 space next to the blood vessel is free (and it's on a different tile)
        blood_vessels = [[0, 4]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

        self.sort_out_halos()

        self.topology.automata[0].grid[0, 3]['contents'] = 'caseum'
        self.topology.automata[0].grid[1, 3]['contents'] = 'caseum'
        self.topology.automata[0].grid[1, 4]['contents'] = 'caseum'
        self.topology.automata[1].grid[1, 0]['contents'] = 'caseum'
        self.topology.automata[0].caseum += [[0, 3], [1, 3], [1, 4]]
        self.topology.automata[1].caseum.append([1,0])

        self.topology.automata[0].update()
        self.topology.automata[1].update()
        self.assertTrue(self.topology.automata[0].potential_events[0].dependant_addresses[0] == [0, 5])

        event = self.topology.automata[0].potential_events[0]

        new_address = self.topology.local_to_local(0, event.dependant_addresses[0], 1)
        self.assertTrue(new_address == [0, 0])

        new_event = event.clone([new_address])

        self.topology.automata[0].process_events([event])
        self.topology.automata[1].process_events([new_event])

        self.assertEqual(len(self.topology.automata[0].t_cells), 0)
        self.assertEqual(len(self.topology.automata[1].t_cells), 1)

    def test_t_cell_recruit_negative_global_bacteria(self):
        self.sort_out_halos()
        # Alter the parameter
        for i in range(4):
            self.topology.automata[i].parameters['bacteria_threshold_for_t_cells'] = 1000000
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_t_cell_recruit_negative_probability(self):
        self.sort_out_halos()
        # Alter the parameter
        for i in range(4):
            self.topology.automata[i].parameters['t_cell_recruitment_probability'] = 0
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_t_cell_recruit_negative_scale(self):
        self.sort_out_halos()
        # Alter the parameter
        for i in range(4):
            self.topology.automata[i].parameters['chemokine_scale_for_t_cell_recruitment'] = 101
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)


class MacrophageRecruitmentTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['time_step'] = 0.001
        params['chemokine_diffusion'] = 0.05
        params['chemokine_decay'] = 0.347
        params['chemotherapy_schedule1_start'] = 100
        params['chemotherapy_schedule2_start'] = 200
        params['t_cell_movement_time'] = 1
        params['bacteria_threshold_for_t_cells'] = 1000
        params['interval_to_record_results'] = 1000

        params['macrophage_recruitment_probability'] = 100
        params['chemokine_scale_for_macrophage_recruitment'] = -1

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        blood_vessels = [[3, 3]]
        fast_bacteria = [[9, 9], [8, 9]]
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        bacteria_global = 0
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
            bacteria_global += len(i.bacteria)
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])
            self.topology.automata[i].set_global_bacteria_number(bacteria_global)

    def test_macrophage_recruitment(self):
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.RecruitMacrophage))

        address = event.dependant_addresses[0]
        distance = math.fabs(address[0] - 3) + math.fabs(address[1] - 3)
        self.assertEqual(distance, 1)

    def test_macrophage_recruit_process(self):
        self.sort_out_halos()
        self.topology.automata[0].update()
        events = self.topology.automata[0].potential_events
        self.topology.automata[0].process_events(events)

        self.assertEqual(len(self.topology.automata[0].macrophages), 1)
        self.assertTrue(isinstance(self.topology.automata[0].get_attribute(events[0].dependant_addresses[0], 'contents'),
                                   TB_Model.Macrophage))

    def test_macrophage_recruited_across_boundary(self):

        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 1.0
        params['chemotherapy_diffusion'] = 0.75
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['caseum_threshold_to_reduce_diffusion'] = 100
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 2.4
        params['time_step'] = 0.001
        params['chemokine_diffusion'] = 0.05
        params['chemokine_decay'] = 0.347
        params['chemotherapy_schedule1_start'] = 100
        params['chemotherapy_schedule2_start'] = 200
        params['t_cell_movement_time'] = 1
        params['bacteria_threshold_for_t_cells'] = 1000
        params['interval_to_record_results'] = 1000

        params['macrophage_recruitment_probability'] = 100
        params['chemokine_scale_for_macrophage_recruitment'] = -1

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        blood_vessels = [[0, 4]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

        self.sort_out_halos()

        self.topology.automata[0].grid[0, 3]['contents'] = 'caseum'
        self.topology.automata[0].grid[1, 3]['contents'] = 'caseum'
        self.topology.automata[0].grid[1, 4]['contents'] = 'caseum'
        self.topology.automata[1].grid[1, 0]['contents'] = 'caseum'
        self.topology.automata[0].caseum += [[0, 3],[1, 3],[1, 4]]
        self.topology.automata[1].caseum.append([1,0])

        self.topology.automata[0].update()
        self.assertTrue(self.topology.automata[0].potential_events[0].dependant_addresses[0] == [0, 5])

        self.topology.automata[1].update()

        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.RecruitMacrophage))

        new_address = self.topology.local_to_local(0, event.dependant_addresses[0], 1)
        self.assertTrue(new_address == [0, 0])

        new_event = event.clone([new_address])

        self.topology.automata[0].process_events([event])
        self.topology.automata[1].process_events([new_event])

        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(len(self.topology.automata[1].macrophages), 1)

    def test_macrophage_recruit_negative_probability(self):
        self.sort_out_halos()
        for i in self.topology.automata:
            i.parameters['macrophage_recruitment_probability'] = 0

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_recruit_negative_scale(self):
        self.sort_out_halos()
        for i in self.topology.automata:
            i.parameters['chemokine_scale_for_macrophage_recruitment'] = 101

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)


class ChemotherapyKillsBacteriaTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 0.0
        params['chemotherapy_diffusion'] = 0.0
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['spatial_step'] = 0.2
        params['oxygen_from_source'] = 0.0
        params['time_step'] = 0.001
        params['chemokine_diffusion'] = 0.0
        params['chemokine_decay'] = 0.0
        params['oxygen_uptake_from_bacteria'] = 0.0
        params['chemokine_from_bacteria'] = 0.0
        params['macrophage_recruitment_probability'] = 0
        params['bacteria_threshold_for_t_cells'] = 1000
        params['bacteria_replication_fast_upper'] = 99999
        params['bacteria_replication_fast_lower'] = 99998
        params['bacteria_replication_slow_upper'] = 99999
        params['bacteria_replication_slow_lower'] = 99998
        params['interval_to_record_results'] = 1000

        # Scale = 0, so any chemo will kill bacteria (will change later for negative tests)
        params['chemotherapy_scale_for_kill_fast_bacteria'] = 0
        params['chemotherapy_scale_for_kill_slow_bacteria'] = 0
        # No chemo - will manually add chemotherapy to necessary squares
        params['chemotherapy_schedule1_start'] = 99
        params['chemotherapy_schedule1_end'] = 100
        params['chemotherapy_schedule2_start'] = 200

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        blood_vessels = [[3, 3]]
        fast_bacteria = [[1, 1]]
        slow_bacteria = [[2, 7]]
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        max_chemo = 0
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
            max_chemo = max(max_chemo, i.max_chemotherapy_local)

        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])
            self.topology.automata[i].set_max_chemotherapy_global(max_chemo)

    def test_chemo_kills_bacteria(self):

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1
        self.topology.automata[1].grid[2, 2]['chemotherapy'] = 1
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2
        self.topology.automata[0].max_chemotherapy_local = 1
        self.topology.automata[1].max_chemotherapy_local = 1

        # So max chemo is 2, and 1,1 / 2,2 have 50%.

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillBacterium))
        self.assertTrue(event.dependant_addresses[0] == [1, 1])

        self.topology.automata[1].update()
        self.assertEqual(len(self.topology.automata[1].potential_events), 1)
        self.assertTrue(isinstance(self.topology.automata[1].potential_events[0], TB_Model.ChemoKillBacterium))
        event = self.topology.automata[1].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillBacterium))
        self.assertTrue(event.dependant_addresses[0] == [2, 2])

    def test_chemokillbacteria_negative_both(self):

        for i in self.topology.automata:
            i.parameters['chemotherapy_scale_for_kill_fast_bacteria'] = 70
            i.parameters['chemotherapy_scale_for_kill_slow_bacteria'] = 70

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1
        self.topology.automata[1].grid[2, 2]['chemotherapy'] = 1
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2
        self.topology.automata[0].max_chemotherapy_local = 1
        self.topology.automata[1].max_chemotherapy_local = 1
        self.topology.automata[3].max_chemotherapy_local = 2

        # So max chemo is 2, and 1,1 / 2,2 have 50%.

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

        self.topology.automata[1].update()
        self.assertEqual(len(self.topology.automata[1].potential_events), 0)

    def test_chemokillsbacteria_fast_but_not_slow(self):
        for i in self.topology.automata:
            i.parameters['chemotherapy_scale_for_kill_fast_bacteria'] = 20.0
            i.parameters['chemotherapy_scale_for_kill_slow_bacteria'] = 70.0

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[1].grid[2, 2]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[1].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillBacterium))
        self.assertTrue(event.dependant_addresses[0] == [1, 1])

        self.topology.automata[1].update()
        self.assertEqual(len(self.topology.automata[1].potential_events), 0)

    def test_chemokillbacteria_process(self):
        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1
        self.topology.automata[1].grid[2, 2]['chemotherapy'] = 1
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2
        self.topology.automata[0].max_chemotherapy_local = 1
        self.topology.automata[1].max_chemotherapy_local = 1

        # So max chemo is 2, and 1,1 / 2,2 have 50%.

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event_0 = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event_0, TB_Model.ChemoKillBacterium))

        self.topology.automata[1].update()
        self.assertEqual(len(self.topology.automata[1].potential_events), 1)
        self.assertTrue(isinstance(self.topology.automata[1].potential_events[0], TB_Model.ChemoKillBacterium))
        event_1 = self.topology.automata[1].potential_events[0]
        self.assertTrue(isinstance(event_1, TB_Model.ChemoKillBacterium))

        # Now process
        self.topology.automata[0].process_events([event_0])
        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

        self.topology.automata[1].process_events([event_1])
        self.assertEqual(len(self.topology.automata[1].bacteria), 0)
        self.assertEqual(self.topology.automata[1].grid[2, 2]['contents'], 0.0)


class ChemotherapyKillsMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        params['initial_oxygen'] = 1.5
        params['oxygen_diffusion'] = 0.0
        params['chemotherapy_diffusion'] = 0.0
        params['caseum_distance_to_reduce_diffusion'] = 2
        params['spatial_step'] = 0.2
        params['time_step'] = 0.001
        params['oxygen_from_source'] = 0.0
        params['chemokine_diffusion'] = 0.0
        params['chemokine_decay'] = 0.0
        params['chemokine_from_macrophage'] = 0.0
        params['bacteria_threshold_for_t_cells'] = 1000
        params['macrophage_recruitment_probability'] = 0
        params['chemokine_scale_for_macrophage_activation'] = 101
        params['resting_macrophage_age_limit'] = 1000000
        params['resting_macrophage_movement_time'] = 1000000
        params['active_macrophage_age_limit'] = 999
        params['active_macrophage_movement_time'] = 999
        params['infected_macrophage_age_limit'] = 999
        params['infected_macrophage_movement_time'] = 999
        params['chronically_infected_macrophage_age_limit'] = 999
        params['chronically_infected_macrophage_movement_time'] = 999
        params['chemokine_scale_for_macrophage_deactivation'] = -1
        params['bacteria_to_turn_chronically_infected'] = 99
        params['bacteria_to_burst_macrophage'] = 99
        params['interval_to_record_results'] = 1000

        params['chemotherapy_scale_for_kill_macrophage'] = 0

        params['chemotherapy_schedule1_start'] = 99
        params['chemotherapy_schedule2_start'] = 200

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        blood_vessels = [[3, 3]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[1, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        max_chemo = 0
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
            max_chemo = max(max_chemo, i.max_chemotherapy_local)

        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])
            self.topology.automata[i].set_max_chemotherapy_global(max_chemo)

    def test_chemokillmacrophage_infected(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'

        # Only 1 macrophage and it's resting

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillMacrophage))
        address = event.dependant_addresses[0]
        self.assertTrue(address == [1, 1])

    def test_chemokillmacrophage_chronic_infected(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'

        # Only 1 macrophage and it's resting

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillMacrophage))
        address = event.dependant_addresses[0]
        self.assertTrue(address == [1, 1])

    def test_process_chemokillmacrophage(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'

        # Only 1 macrophage and it's resting

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.ChemoKillMacrophage))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

    def test_chemokillmacrophage_negative_resting(self):

        # Only 1 macrophage and it's resting

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_chemokillmacrophage_negative_active(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        # Only 1 macrophage and it's active

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_chemokillmacrophage_negative_scale(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'

        for i in self.topology.automata:
            i.parameters['chemotherapy_scale_for_kill_macrophage'] = 75.

        # Only 1 macrophage and it's resting

        self.topology.automata[0].grid[1, 1]['chemotherapy'] = 1.0
        self.topology.automata[3].grid[4, 4]['chemotherapy'] = 2.0
        self.topology.automata[0].max_chemotherapy_local = 1.0
        self.topology.automata[3].max_chemotherapy_local = 2.0

        # So 1,1 has scale of 50%

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)


class TCellDeathTestCase(unittest.TestCase):
    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['macrophage_recruitment_probability'] = 0
        params['chemotherapy_scale_for_kill_macrophage'] = 0
        params['t_cell_random_move_probability'] = 100
        params['interval_to_record_results'] = 1000

        params['bacteria_threshold_for_t_cells'] = 0
        params['t_cell_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['t_cell_age_threshold'] = 2
        params['time_step'] = 2

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[3, 3]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_t_cell_death(self):

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellDeath))
        address = event.dependant_addresses[0]
        self.assertTrue(address == [1, 1])

    def test_process_tcell_death(self):

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellDeath))

        # PROCESS
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].t_cells), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

    def test_tcell_death_negative_bacteria_threshold(self):
        # Don't think this makes sense, but based on Ruth's code (i.e. why does number of bacteria stop/cause
        # t-cell death)

        for i in self.topology.automata:
            i.parameters['bacteria_threshold_for_t_cells'] = 1

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_tcell_death_negative_movement_time(self):
        # Don't think this makes sense, but based on Ruth's code (i.e. why does movement time for t-cells stop/cause
        # t-cell death)

        for i in self.topology.automata:
            i.parameters['t_cell_movement_time'] = 15

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_tcell_death_negative_age_threshold(self):
        for i in self.topology.automata:
            i.parameters['t_cell_age_threshold'] = 999

        # Set seed - forces random age to be 320
        np.random.seed(10)

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        # Cells either move or die when time % movement_time = 0 so check it's a move
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellMovement))

        # Just in case
        np.random.seed(None)


class TCellMovementTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['macrophage_recruitment_probability'] = 0
        params['chemotherapy_scale_for_kill_macrophage'] = 0
        params['interval_to_record_results'] = 1000

        params['bacteria_threshold_for_t_cells'] = 0
        params['t_cell_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['t_cell_age_threshold'] = 1000
        params['time_step'] = 2
        params['t_cell_random_move_probability'] = 100

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[3, 3]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_t_cell_movement_random(self):

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Forces [0,1] as 'random' choice
        np.random.seed(10)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellMovement))
        self.assertTrue(event.dependant_addresses[0] == [1, 1])
        self.assertSequenceEqual(event.dependant_addresses[1], [0, 1])

        # Just in case
        np.random.seed(None)

    def test_t_cell_movement_max_chemokine(self):

        # Never random
        for i in self.topology.automata:
            i.parameters['t_cell_random_move_probability'] = 0.0

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[0, 0]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellMovement))

        self.assertTrue(event.dependant_addresses[0] == [1, 1])
        self.assertSequenceEqual(event.dependant_addresses[1], [0, 0])

    def test_t_cell_movement_process(self):
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Forces [0,1] as 'random' choice
        np.random.seed(10)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellMovement))
        self.assertTrue(event.dependant_addresses[0] == [1, 1])
        self.assertSequenceEqual(event.dependant_addresses[1], [0, 1])

        # Now process
        self.topology.automata[0].process_events([event])

        self.assertTrue(isinstance(self.topology.automata[0].grid[0, 1]['contents'], TB_Model.TCell))
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

    def test_t_cell_movement_across_boundary_process(self):

        # Never random
        for i in self.topology.automata:
            i.parameters['t_cell_random_move_probability'] = 0.0

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([4, 4])
        self.topology.automata[0].grid[4, 4]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[1].grid[4, 0]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)
        self.topology.automata[1].set_max_chemokine_global(99.9)

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.topology.automata[1].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellMovement))

        self.assertSequenceEqual(event.dependant_addresses[0], [4, 4])
        self.assertSequenceEqual(event.dependant_addresses[1], [4, 5])

        # Now process
        new_addresses = [self.topology.local_to_local(0, event.dependant_addresses[0], 1),
                         self.topology.local_to_local(0, event.dependant_addresses[1], 1)]

        self.assertSequenceEqual(new_addresses[0], [4, -1])
        self.assertSequenceEqual(new_addresses[1], [4, 0])

        new_event = event.clone(new_addresses)
        self.assertSequenceEqual(new_event.dependant_addresses[0], [4, -1])
        self.assertSequenceEqual(new_event.dependant_addresses[1], [4, 0])

        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[4, 4]['contents'], 0.0)

        self.topology.automata[1].process_events([new_event])
        self.assertTrue(isinstance(self.topology.automata[1].grid[4, 0]['contents'], TB_Model.TCell))

    def test_t_cell_movement_negative_number_bacteria(self):

        for i in self.topology.automata:
            i.parameters['bacteria_threshold_for_t_cells'] = 1

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_t_cell_movement_negative_movement_time(self):

        for i in self.topology.automata:
            i.parameters['t_cell_movement_time'] = 99

        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)


class TCellKillsMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['macrophage_recruitment_probability'] = 0
        params['chemotherapy_scale_for_kill_macrophage'] = 0
        params['chemokine_from_macrophage'] = 0
        params['resting_macrophage_age_limit'] = 999
        params['resting_macrophage_movement_time'] = 999
        params['active_macrophage_age_limit'] = 999
        params['active_macrophage_movement_time'] = 999
        params['infected_macrophage_age_limit'] = 999
        params['infected_macrophage_movement_time'] = 999
        params['chronically_infected_macrophage_age_limit'] = 999
        params['chronically_infected_macrophage_movement_time'] = 999
        params['bacteria_to_turn_chronically_infected'] = 99
        params['chemokine_scale_for_macrophage_deactivation'] = -1
        params['chemokine_scale_for_macrophage_activation'] = 101
        params['bacteria_to_burst_macrophage'] = 99
        params['interval_to_record_results'] = 1000

        params['bacteria_threshold_for_t_cells'] = 0
        params['t_cell_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['t_cell_age_threshold'] = 1000
        params['time_step'] = 2
        params['t_cell_random_move_probability'] = 0
        params['t_cell_kills_macrophage_probability'] = 100

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[1, 2]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_tcell_kill_macrophage_infected_event(self):

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellKillsMacrophage))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])

    def test_tcell_kill_macrophage_chronic_infected_event(self):
        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'chronically_infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellKillsMacrophage))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])

    def test_tcell_kill_macrophage_process(self):
        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellKillsMacrophage))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertEqual(self.topology.automata[0].grid[1, 2]['contents'], 'caseum')
        self.assertTrue([1,2] in self.topology.automata[0].caseum)

    def test_tcell_kill_macrophage_negative_resting(self):

        np.random.seed(100)

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'resting'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_tcell_kill_macrophage_negative_active(self):
        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'active'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_tcell_kill_macrophage_negative_probability(self):

        np.random.seed(100)

        for i in self.topology.automata:
            i.parameters['t_cell_kills_macrophage_probability'] = 0

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([1, 1])
        self.topology.automata[0].grid[1, 1]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)

        self.topology.automata[0].macrophages[0].state = 'infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_tcell_kills_macrophage_across_boundary_process(self):
        params = dict()
        params['max_depth'] = 3
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
        params['macrophage_recruitment_probability'] = 0
        params['chemotherapy_scale_for_kill_macrophage'] = 0
        params['chemokine_from_macrophage'] = 0
        params['resting_macrophage_age_limit'] = 999
        params['resting_macrophage_movement_time'] = 999
        params['active_macrophage_age_limit'] = 999
        params['active_macrophage_movement_time'] = 999
        params['infected_macrophage_age_limit'] = 999
        params['infected_macrophage_movement_time'] = 999
        params['chronically_infected_macrophage_age_limit'] = 999
        params['chronically_infected_macrophage_movement_time'] = 999
        params['interval_to_record_results'] = 1000
        params['bacteria_to_turn_chronically_infected'] = 99

        params['bacteria_threshold_for_t_cells'] = 0
        params['t_cell_recruitment_probability'] = 0
        params['t_cell_movement_time'] = 1
        params['t_cell_age_threshold'] = 1000
        params['time_step'] = 2
        params['t_cell_random_move_probability'] = 0
        params['t_cell_kills_macrophage_probability'] = 100

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 5]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

        # Add a t-cell to [1,1]
        t_cell = TB_Model.TCell([4, 4])
        self.topology.automata[0].grid[4, 4]['contents'] = t_cell
        self.topology.automata[0].t_cells.append(t_cell)

        # Make [0,0] have the highest chemokine
        self.topology.automata[1].grid[4, 0]['chemokine'] = 99.9
        self.topology.automata[0].set_max_chemokine_global(99.9)
        self.topology.automata[1].set_max_chemokine_global(99.9)

        self.topology.automata[1].macrophages[0].state = 'infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.TCellKillsMacrophage))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

        self.topology.automata[1].update()

        new_addresses = [self.topology.local_to_local(0, addresses[0], 1),
                         self.topology.local_to_local(0, addresses[1], 1)]

        self.assertSequenceEqual(new_addresses[0], [4, -1])
        self.assertSequenceEqual(new_addresses[1], [4, 0])

        new_event = event.clone(new_addresses)
        self.assertSequenceEqual(new_event.dependant_addresses[0], [4, -1])
        self.assertSequenceEqual(new_event.dependant_addresses[1], [4, 0])

        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[4, 4]['contents'], 0.0)
        self.assertEqual(len(self.topology.automata[0].t_cells), 0)

        self.topology.automata[1].process_events([new_event])
        self.assertEqual(self.topology.automata[1].grid[4, 0]['contents'], 'caseum')
        self.assertTrue([4,0] in self.topology.automata[1].caseum)
        self.assertEqual(len(self.topology.automata[1].macrophages), 0)


class MacrophageDeathTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['chemokine_scale_for_macrophage_activation'] = 101
        params['chemotherapy_scale_for_kill_macrophage'] = 101
        params['resting_macrophage_movement_time'] = 999
        params['active_macrophage_movement_time'] = 999
        params['infected_macrophage_movement_time'] = 999
        params['chronically_infected_macrophage_movement_time'] = 999
        params['chemokine_scale_for_macrophage_deactivation'] = -1
        params['bacteria_to_turn_chronically_infected'] = 99
        params['bacteria_to_burst_macrophage'] = 99
        params['interval_to_record_results'] = 1000

        params['time_step'] = 1
        params['resting_macrophage_age_limit'] = 2
        params['active_macrophage_age_limit'] = 4
        params['infected_macrophage_age_limit'] = 2
        params['chronically_infected_macrophage_age_limit'] = 2

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[1, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_macrophage_death_resting(self):

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [1, 1])

    def test_macrophage_death_active(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[1, 1]['contents'].age = 999

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [1, 1])

    def test_macrophage_death_infected(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [1, 1])

    def test_macrophage_death_chronically_infected(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [1, 1])

    def test_macrophage_resting_death_negative_age(self):

        self.topology.automata[0].parameters['resting_macrophage_age_limit'] = 200

        # Seed forces age to 9
        np.random.seed(100)

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_active_death_negative_age(self):

        self.topology.automata[0].parameters['active_macrophage_age_limit'] = 200

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[1, 1]['contents'].age = 2

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_infected_death_negative_age(self):

        self.topology.automata[0].parameters['infected_macrophage_age_limit'] = 200

        # Seed forces age to 9
        np.random.seed(100)

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_chronically_infected_death_negative_age(self):

        self.topology.automata[0].parameters['chronically_infected_macrophage_age_limit'] = 200

        # Seed forces age to 9
        np.random.seed(100)

        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_death_resting_process(self):
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

    def test_macrophage_death_active_process(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[1, 1]['contents'].age = 999

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)

    def test_macrophage_death_infected_process(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 'caseum')
        self.assertTrue([1,1] in self.topology.automata[0].caseum)

    def test_macrophage_death_chronically_infected_process(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageDeath))

        # Now process
        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 'caseum')
        self.assertTrue([1, 1] in self.topology.automata[0].caseum)


class MacrophageMovementTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['chemokine_scale_for_macrophage_activation'] = 101
        params['chemotherapy_scale_for_kill_macrophage'] = 101
        params['chemokine_scale_for_macrophage_deactivation'] = -1
        params['bacteria_to_turn_chronically_infected'] = 99
        params['bacteria_to_burst_macrophage'] = 99
        params['interval_to_record_results'] = 1000

        params['time_step'] = 1
        params['resting_macrophage_age_limit'] = 999
        params['active_macrophage_age_limit'] = 999
        params['infected_macrophage_age_limit'] = 999
        params['chronically_infected_macrophage_age_limit'] = 999
        params['resting_macrophage_movement_time'] = 1
        params['active_macrophage_movement_time'] = 1
        params['infected_macrophage_movement_time'] = 1
        params['chronically_infected_macrophage_movement_time'] = 1
        params['prob_resting_macrophage_random_move'] = 100
        params['minimum_chemokine_for_resting_macrophage_movement'] = 101

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[1, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_macrophage_resting_movement_random_through_prob(self):
        self.sort_out_halos()

        # Random choice is [1,0]
        np.random.seed(100)

        # Set [1,2] as max chemokine cell - make sure it doesn't go here though (cause it's a random move)
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99
        self.topology.automata[0].max_chemokine_global = 99

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageMovement)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 0])

    def test_macrophage_resting_movement_random_through_lack_of_chemokine(self):

        # Random choice is [1,0]
        np.random.seed(100)

        # Not random (sort of)
        self.topology.automata[0].parameters['prob_resting_macrophage_random_move'] = 0

        # Set [1,2] as max chemokine cell - make sure it doesn't go here though (cause it's scale is too low)
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99999.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageMovement)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 0])

    def test_macrophage_resting_movement_max_chemokine(self):
        # Not random
        self.topology.automata[0].parameters['prob_resting_macrophage_random_move'] = 0
        # Any chemokine forces a move there
        self.topology.automata[0].parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0

        # Set [1,0] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 0]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageMovement)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 0])

    def test_macrophage_resting_movement_negative_movement_time(self):
        self.topology.automata[0].parameters['resting_macrophage_movement_time'] = 99

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_resting_movement_across_boundary(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 4]]

        # Not random
        self.parameters['prob_resting_macrophage_random_move'] = 0
        # Any chemokine forces a move there
        self.parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        self.assertFalse(event.internal)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

    def test_macrophage_active_movement_max_chemokine(self):
        # Active always moves to max chemokine
        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[2, 1]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [2, 1])

    def test_macrophage_active_movement_across_boundary(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 4]]

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100
        self.topology.automata[0].grid[4, 4]['contents'].state = 'active'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        self.assertFalse(event.internal)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

    def test_macrophage_active_movement_negative_movement_time(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].parameters['active_macrophage_movement_time'] = 99

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_infected_movement_max_chemokine(self):
        # Active always moves to max chemokine
        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.topology.automata[0].grid[2, 1]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [2, 1])

    def test_macrophage_infected_movement_across_boundary(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 4]]

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100
        self.topology.automata[0].grid[4, 4]['contents'].state = 'infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        self.assertFalse(event.internal)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

    def test_macrophage_infected_movement_negative_movement_time(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.topology.automata[0].parameters['infected_macrophage_movement_time'] = 99

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_chronically_infected_movement_max_chemokine(self):
        # Active always moves to max chemokine
        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[2, 1]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [2, 1])

    def test_macrophage_chronically_infected_movement_across_boundary(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 4]]

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100
        self.topology.automata[0].grid[4, 4]['contents'].state = 'chronically_infected'

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        self.assertFalse(event.internal)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

    def test_macrophage_chronically_infected_movement_negative_movement_time(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.topology.automata[0].parameters['chronically_infected_macrophage_movement_time'] = 99

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_process_macrophage_movement(self):
        # Not random
        self.topology.automata[0].parameters['prob_resting_macrophage_random_move'] = 0
        # Any chemokine forces a move there
        self.topology.automata[0].parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0

        # Set [1,0] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 0]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageMovement)

        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 0]['contents'], TB_Model.Macrophage))

    def test_process_macrophage_movement_across_boundaries(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[4, 4]]

        # Not random
        self.parameters['prob_resting_macrophage_random_move'] = 0
        # Any chemokine forces a move there
        self.parameters['minimum_chemokine_for_resting_macrophage_movement'] = 0

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageMovement))
        self.assertFalse(event.internal)

        self.topology.automata[1].update()

        new_addresses = [self.topology.local_to_local(0, event.dependant_addresses[0], 1),
                         self.topology.local_to_local(0, event.dependant_addresses[1], 1)]

        new_event = event.clone(new_addresses)

        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[4, 4]['contents'], 0.0)

        self.topology.automata[1].process_events([new_event])
        self.assertEqual(len(self.topology.automata[1].macrophages), 1)
        self.assertTrue(isinstance(self.topology.automata[1].grid[4, 0]['contents'], TB_Model.Macrophage))


class MacrophageKillsBacteria(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['chemokine_scale_for_macrophage_activation'] = 101
        params['chemotherapy_scale_for_kill_macrophage'] = 101
        params['oxygen_uptake_from_bacteria'] = 0
        params['chemokine_from_bacteria'] = 0
        params['bacteria_replication_fast_upper'] = 999
        params['bacteria_replication_fast_lower'] = 998
        params['bacteria_replication_slow_upper'] = 999
        params['bacteria_replication_slow_lower'] = 998
        params['chemotherapy_scale_for_kill_fast_bacteria'] = 101
        params['chemotherapy_scale_for_kill_slow_bacteria'] = 101
        params['chemokine_scale_for_macrophage_deactivation'] = -1
        params['bacteria_to_turn_chronically_infected'] = 99
        params['bacteria_to_burst_macrophage'] = 99
        params['interval_to_record_results'] = 1000

        params['time_step'] = 1
        params['resting_macrophage_age_limit'] = 999
        params['active_macrophage_age_limit'] = 999
        params['infected_macrophage_age_limit'] = 999
        params['chronically_infected_macrophage_age_limit'] = 999
        params['resting_macrophage_movement_time'] = 1
        params['active_macrophage_movement_time'] = 1
        params['infected_macrophage_movement_time'] = 1
        params['chronically_infected_macrophage_movement_time'] = 1
        params['prob_resting_macrophage_random_move'] = 0
        params['minimum_chemokine_for_resting_macrophage_movement'] = 0

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[8, 8]]
        fast_bacteria = [[1, 2]]
        slow_bacteria = []
        macrophages = [[1, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_resting_macrophage_kill_bacteria(self):
        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])
        self.assertTrue(event.internal)

    def test_resting_macrophage_kill_bacteria_across_boundary(self):

        blood_vessels = [[8, 8]]
        fast_bacteria = [[4, 5]]
        slow_bacteria = []
        macrophages = [[4, 4]]

        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[1].grid[4, 0]['chemokine'] = 100
        self.topology.automata[0].max_chemokine_global = 100
        self.topology.automata[1].max_chemokine_global = 100

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageKillsBacterium))
        self.assertFalse(event.internal)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])

    def test_active_macrophage_kill_fast_bacteria(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_fast_bacteria'] = 100

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])
        self.assertTrue(event.internal)

    def test_active_macrophage_kill_slow_bacteria(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[1, 2]['contents'].metabolism = 'slow'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_slow_bacteria'] = 100

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])
        self.assertTrue(event.internal)

    def test_active_macrophage_kill_fast_bacteria_negative_probability(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_fast_bacteria'] = 0

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_active_macrophage_kill_slow_bacteria_negative_probability(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[1, 2]['contents'].metabolism = 'slow'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_slow_bacteria'] = 0

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_active_macrophage_kill_bacteria_across_boundary(self):

        blood_vessels = [[8, 8]]
        fast_bacteria = [[4, 5]]
        slow_bacteria = []
        macrophages = [[4, 4]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[0].grid[4, 4]['contents'].state = 'active'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_fast_bacteria'] = 100

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[1].grid[4, 0]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9
        self.topology.automata[1].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [4, 4])
        self.assertSequenceEqual(addresses[1], [4, 5])
        self.assertFalse(event.internal)

    def test_infected_macrophage_kill_bacteria(self):
        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])
        self.assertTrue(event.internal)

    def test_chronically_infected_macrophage_kill_bacteria(self):
        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)
        addresses = event.dependant_addresses
        self.assertSequenceEqual(addresses[0], [1, 1])
        self.assertSequenceEqual(addresses[1], [1, 2])
        self.assertTrue(event.internal)

    def test_process_resting_macrophage_kills_bacteria(self):
        # Macrophage kills, it's internal bactera becomes 1

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(len(self.topology.automata[0].macrophages), 1)

        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 2]['contents'], TB_Model.Macrophage))
        self.assertEqual(self.topology.automata[0].grid[1, 2]['contents'].intracellular_bacteria, 1)

    def test_process_active_macrophage_kills_bacteria(self):
        # Macrophage kills, it's internal bactera becomes 1
        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].parameters['prob_active_macrophage_kill_fast_bacteria'] = 100

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(len(self.topology.automata[0].macrophages), 1)

        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 2]['contents'], TB_Model.Macrophage))
        self.assertEqual(self.topology.automata[0].grid[1, 2]['contents'].intracellular_bacteria, 0)

    def test_process_infected_macrophage_kills_bacteria(self):
        # Macrophage kills, it's internal bactera becomes 1
        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.topology.automata[0].grid[1, 1]['contents'].intracellular_bacteria = 4

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(len(self.topology.automata[0].macrophages), 1)

        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 2]['contents'], TB_Model.Macrophage))
        self.assertEqual(self.topology.automata[0].grid[1, 2]['contents'].intracellular_bacteria, 5)

    def test_process_chronically_infected_macrophage_kills_bacteria(self):
        # Macrophage kills, it's internal bactera becomes 1
        self.topology.automata[0].grid[1, 1]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[1, 1]['contents'].intracellular_bacteria = 14

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[0].grid[1, 2]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.assertEqual(len(self.topology.automata[0].bacteria), 1)

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].bacteria), 0)
        self.assertEqual(len(self.topology.automata[0].macrophages), 1)

        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'], 0.0)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 2]['contents'], TB_Model.Macrophage))
        self.assertEqual(self.topology.automata[0].grid[1, 2]['contents'].intracellular_bacteria, 15)

    def test_process_macrophage_kill_bacteria_across_boundary(self):

        blood_vessels = [[8, 8]]
        fast_bacteria = [[4, 5]]
        slow_bacteria = []
        macrophages = [[4, 4]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        # Set [1,2] as max chemokine cell - make sure it goes here
        self.topology.automata[1].grid[4, 0]['chemokine'] = 99.9
        self.topology.automata[0].max_chemokine_global = 99.9
        self.topology.automata[1].max_chemokine_global = 99.9

        self.sort_out_halos()

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(event, TB_Model.MacrophageKillsBacterium)

        self.topology.automata[1].update()

        new_addresses = [self.topology.local_to_local(0, event.dependant_addresses[0], 1),
                         self.topology.local_to_local(0, event.dependant_addresses[1], 1)]

        new_event = event.clone(new_addresses)

        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[4, 4]['contents'], 0.0)
        self.assertEqual(len(self.topology.automata[0].macrophages), 0)

        self.topology.automata[1].process_events([new_event])
        self.assertTrue(isinstance(self.topology.automata[1].grid[4, 0]['contents'], TB_Model.Macrophage))
        self.assertEqual(len(self.topology.automata[1].macrophages), 1)


class MacrophageChangesState(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['interval_to_record_results'] = 1000

        params['chemokine_scale_for_macrophage_activation'] = 0
        params['chemokine_scale_for_macrophage_deactivation'] = 100
        params['time_step'] = 0.5

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[1, 1]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_resting_to_active(self):
        self.topology.automata[0].grid[1, 1]['chemokine'] = 100.0
        self.topology.automata[0].max_chemokine_global = 100.0

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageChangesState))
        self.assertEqual(event.new_state, 'active')

    def test_resting_to_active_negative_scale(self):
        self.topology.automata[0].grid[1, 1]['chemokine'] = 40.0
        self.topology.automata[0].max_chemokine_global = 100.0
        self.topology.automata[0].parameters['chemokine_scale_for_macrophage_activation'] = 50

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_resting_to_infected(self):
        self.topology.automata[0].grid[1, 1]['chemokine'] = 70.0
        self.topology.automata[0].max_chemokine_global = 100.0
        self.topology.automata[0].parameters['chemokine_scale_for_macrophage_activation'] = 50
        self.topology.automata[0].grid[1, 1]['contents'].intracellular_bacteria = 1

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageChangesState))
        self.assertEqual(event.new_state, 'infected')

    def test_active_to_resting(self):

        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'

        self.topology.automata[0].grid[1, 1]['chemokine'] = 20.0
        self.topology.automata[0].max_chemokine_global = 100.0

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageChangesState))
        self.assertEqual(event.new_state, 'resting')

    def test_active_to_resting_negative_scale(self):
        self.topology.automata[0].grid[1, 1]['contents'].state = 'active'
        self.topology.automata[0].parameters['chemokine_scale_for_macrophage_deactivation'] = 50

        self.topology.automata[0].grid[1, 1]['chemokine'] = 50.1
        self.topology.automata[0].max_chemokine_global = 100.0

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_infected_to_chronically(self):
        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.topology.automata[0].grid[1, 1]['contents'].intracellular_bacteria = 11

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageChangesState))
        self.assertEqual(event.new_state, 'chronically_infected')

    def test_infected_to_chronically_negative_not_enough_bacteria(self):
        self.topology.automata[0].grid[1, 1]['contents'].state = 'infected'
        self.topology.automata[0].grid[1, 1]['contents'].intracellular_bacteria = 4

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_changes_state_process(self):
        self.topology.automata[0].grid[1, 1]['chemokine'] = 100.0
        self.topology.automata[0].max_chemokine_global = 100.0

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageChangesState))

        self.topology.automata[0].process_events([event])
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'].state, 'active')


class BacteriaStateChangeTestCase(unittest.TestCase):
    def setUp(self):
        params = dict()
        params['max_depth'] = 3
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
        params['caseum_threshold_to_reduce_diffusion'] = 100
        params['interval_to_record_results'] = 1000

        params['time_step'] = 4
        params['oxygen_scale_for_metabolism_change_to_slow'] = 100
        params['oxygen_scale_for_metabolism_change_to_fast'] = 0

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[6, 2]]
        fast_bacteria = [[1, 1]]
        slow_bacteria = [[7, 7]]
        macrophages = []
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_fast_to_slow(self):
        self.sort_out_halos()
        self.topology.automata[0].grid[1, 1]['oxygen'] = 1
        self.topology.automata[0].max_oxygen_global = 100

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))
        self.assertEqual(event.type_of_change, 'metabolism')
        self.assertEqual(event.new_value, 'slow')

    def test_fast_to_slow_negative_scale(self):
        self.sort_out_halos()
        self.topology.automata[0].parameters['oxygen_scale_for_metabolism_change_to_slow'] = 50.0

        self.topology.automata[0].grid[1, 1]['oxygen'] = 75.0
        self.topology.automata[0].max_oxygen_global = 100.0

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_slow_to_fast(self):
        self.sort_out_halos()
        self.topology.automata[3].grid[2, 2]['oxygen'] = 1.0
        self.topology.automata[3].max_oxygen_global = 100.0

        self.topology.automata[3].update()
        self.assertEqual(len(self.topology.automata[3].potential_events), 1)
        event = self.topology.automata[3].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))
        self.assertEqual(event.type_of_change, 'metabolism')
        self.assertEqual(event.new_value, 'fast')

    def test_slow_to_fast_negative_scale(self):
        self.sort_out_halos()
        self.topology.automata[3].parameters['oxygen_scale_for_metabolism_change_to_fast'] = 50.0

        self.topology.automata[3].grid[2, 2]['oxygen'] = 25.0
        self.topology.automata[3].max_oxygen_global = 100.0

        self.topology.automata[3].update()
        self.assertEqual(len(self.topology.automata[3].potential_events), 0)

    def test_process_bacterium_changes_metabolism(self):
        self.sort_out_halos()
        self.topology.automata[0].grid[1, 1]['oxygen'] = 1
        self.topology.automata[0].max_oxygen_global = 100

        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))

        self.topology.automata[0].process_events([event])
        self.assertEqual(len(self.topology.automata[0].bacteria), 1)
        self.assertTrue(isinstance(self.topology.automata[0].grid[1, 1]['contents'], TB_Model.Bacterium))
        self.assertEqual(self.topology.automata[0].grid[1, 1]['contents'].metabolism, 'slow')

    def test_non_resting_to_resting(self):

        self.topology.automata[0].parameters['bacteria_replication_fast_upper'] = 2
        self.topology.automata[0].parameters['bacteria_replication_fast_lower'] = 1
        self.topology.automata[0].parameters['time_step'] = 1

        # Fill automata 0 with caseum, so forces the bacteria in [1,1] to rest
        for x in range(5):
            for y in range(5):
                if x != 1 or y != 1:
                    self.topology.automata[0].grid[x, y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x,y])

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))
        self.assertEqual(event.type_of_change, 'resting')
        self.assertEqual(event.new_value, True)

    def test_resting_to_non_resting(self):

        self.topology.automata[0].parameters['oxygen_scale_for_metabolism_change_to_slow'] = -1.0

        self.topology.automata[0].grid[1, 1]['contents'].resting = True
        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))
        self.assertEqual(event.type_of_change, 'resting')
        self.assertEqual(event.new_value, False)

    def test_non_resting_to_resting_negative_room(self):
        self.topology.automata[0].parameters['bacteria_replication_fast_upper'] = 2
        self.topology.automata[0].parameters['bacteria_replication_fast_lower'] = 1
        self.topology.automata[0].parameters['time_step'] = 1

        # Fill automata 0 with caseum, so forces the bacteria in [1,1] to rest
        for x in range(5):
            for y in range(5):
                if x != 1 or y != 1:
                    self.topology.automata[0].grid[x, y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x,y])

        self.topology.automata[0].grid[0, 0]['contents'] = 0.0

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        self.assertTrue(isinstance(self.topology.automata[0].potential_events[0], TB_Model.BacteriumReplication))

    def test_resting_to_non_resting_negative_no_room(self):

        self.topology.automata[0].parameters['oxygen_scale_for_metabolism_change_to_slow'] = -1.0

        self.topology.automata[0].grid[1, 1]['contents'].resting = True
        # Fill automata 0 with caseum, so forces the bacteria in [1,1] to rest
        for x in range(5):
            for y in range(5):
                if x != 1 or y != 1:
                    self.topology.automata[0].grid[x, y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x, y])

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_non_resting_to_resting_process(self):

        self.topology.automata[0].parameters['bacteria_replication_fast_upper'] = 2
        self.topology.automata[0].parameters['bacteria_replication_fast_lower'] = 1
        self.topology.automata[0].parameters['time_step'] = 1

        # Fill automata 0 with caseum, so forces the bacteria in [1,1] to rest
        for x in range(5):
            for y in range(5):
                if x != 1 or y != 1:
                    self.topology.automata[0].grid[x, y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x, y])

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))

        self.topology.automata[0].process_events([event])

        self.assertTrue(self.topology.automata[0].grid[1, 1]['contents'].resting)

    def test_resting_to_non_resting_process(self):

        self.topology.automata[0].parameters['oxygen_scale_for_metabolism_change_to_slow'] = -1.0

        self.topology.automata[0].grid[1, 1]['contents'].resting = True
        self.sort_out_halos()
        self.topology.automata[0].update()
        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.BacteriumStateChange))

        self.topology.automata[0].process_events([event])

        self.assertFalse(self.topology.automata[0].grid[1, 1]['contents'].resting)


class MacrophageBurstingTestCase(unittest.TestCase):
    def setUp(self):
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
        params['interval_to_record_results'] = 1000

        params['bacteria_to_burst_macrophage'] = 20

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[6, 2]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[2, 2]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)
        np.random.seed(100)

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_macrophage_bursting(self):

        self.topology.automata[0].grid[2, 2]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[2, 2]['contents'].intracellular_bacteria = 20

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageBursting))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [2, 2])
        self.assertEqual(len(event.bacteria_addresses), 20)

    def test_macrophage_bursting_bacteria_addresses(self):
        self.topology.automata[0].grid[2, 2]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[2, 2]['contents'].intracellular_bacteria = 20

        # A 2x2 neighbourhood has 24 cells, so put caseum in 4 cells - rest will be filled by bacteria
        self.topology.automata[0].grid[0, 0]['contents'] = 'caseum'
        self.topology.automata[0].grid[0, 4]['contents'] = 'caseum'
        self.topology.automata[0].grid[4, 0]['contents'] = 'caseum'
        self.topology.automata[0].grid[4, 4]['contents'] = 'caseum'
        self.topology.automata[0].caseum += [[0,0],[0,4],[4,0],[4,4]]

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageBursting))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [2, 2])
        self.assertEqual(len(event.bacteria_addresses), 20)

        # Check that bacteria have been put into all cells in a 2x2 neighbourhood except the corners (which have caseum)
        self.assertItemsEqual(event.bacteria_addresses,
                              [[0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [2, 0], [2, 1], [2, 3],
                               [2, 4], [3, 0], [3, 1], [3, 2], [3, 3], [3, 4], [4, 1], [4, 2], [4, 3]])

    def test_macrophage_bursting_bacteria_addresses_not_enough_room(self):
        blood_vessels = [[8, 8]]
        fast_bacteria = []
        slow_bacteria = []
        macrophages = [[3, 3]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], self.attributes, self.parameters,
                                                        blood_vessels, fast_bacteria, slow_bacteria, macrophages)

        self.topology.automata[0].grid[3, 3]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[3, 3]['contents'].intracellular_bacteria = 20

        # Fill the grid with so much caseum that there's not enough room for 20 bacteria in a 3x3 neighbourhood
        # Caseum in automaton 0
        for x in range(0, 5):
            for y in range(0, 5):
                if not(x == 3 and y == 3):
                    self.topology.automata[0].grid[x, y]['contents'] = 'caseum'
                    self.topology.automata[0].caseum.append([x,y])
        # Caseum in automaton 1
        for x in range(0, 5):
            for y in range(0, 2):
                self.topology.automata[1].grid[x, y]['contents'] = 'caseum'
                self.topology.automata[1].caseum.append([x, y])

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageBursting))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [3, 3])

        # Check that bacteria have been put into all available cells
        self.assertEqual(len(event.bacteria_addresses), 14)
        self.assertItemsEqual(event.bacteria_addresses,
                              [[5, 0], [5, 1], [5, 2], [5, 3], [5, 4], [5, 5], [5, 6], [6, 0], [6, 1], [6, 2], [6, 3],
                               [6, 4], [6, 5], [6, 6]])

    def test_macrophage_bursting_negative_number_bacteria(self):

        self.topology.automata[0].grid[2, 2]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[2, 2]['contents'].intracellular_bacteria = 19

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 0)

    def test_macrophage_bursting_process(self):

        self.topology.automata[0].grid[2, 2]['contents'].state = 'chronically_infected'
        self.topology.automata[0].grid[2, 2]['contents'].intracellular_bacteria = 20

        # A 2x2 neighbourhood has 24 cells, so put caseum in 4 cells - rest will be filled by bacteria
        self.topology.automata[0].grid[0, 0]['contents'] = 'caseum'
        self.topology.automata[0].grid[0, 4]['contents'] = 'caseum'
        self.topology.automata[0].grid[4, 0]['contents'] = 'caseum'
        self.topology.automata[0].grid[4, 4]['contents'] = 'caseum'
        self.topology.automata[0].caseum += [[0,0], [0,4], [4,0], [4,4]]

        self.sort_out_halos()
        self.topology.automata[0].update()

        self.assertEqual(len(self.topology.automata[0].potential_events), 1)
        event = self.topology.automata[0].potential_events[0]
        self.assertTrue(isinstance(event, TB_Model.MacrophageBursting))
        address = event.dependant_addresses[0]
        self.assertSequenceEqual(address, [2, 2])
        self.assertEqual(len(event.bacteria_addresses), 20)

        # Check that bacteria have been put into all cells in a 2x2 neighbourhood except the corners (which have caseum)
        expected_bacteria_addresses = [[0, 1], [0, 2], [0, 3], [1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [2, 0], [2, 1],
                                       [2, 3], [2, 4], [3, 0], [3, 1], [3, 2], [3, 3], [3, 4], [4, 1], [4, 2], [4, 3]]
        self.assertItemsEqual(event.bacteria_addresses, expected_bacteria_addresses)

        event.impacted_addresses_allowed = event.impacted_addresses_potential

        self.topology.automata[0].process_events([event])

        self.assertEqual(len(self.topology.automata[0].macrophages), 0)
        self.assertEqual(self.topology.automata[0].grid[2, 2]['contents'], 'caseum')
        self.assertTrue([2,2] in self.topology.automata[0].caseum)

        for i in expected_bacteria_addresses:
            cell_contents = self.topology.automata[0].grid[i[0], i[1]]['contents']
            self.assertTrue(isinstance(cell_contents, TB_Model.Bacterium))
            self.assertEqual(cell_contents.metabolism, 'slow')

class OutputWritingTestCase(unittest.TestCase):

    def setUp(self):
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
        params['interval_to_record_results'] = 1000
        params['bacteria_to_burst_macrophage'] = 20
        params['macrophage_recruitment_probability'] = 0

        atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
                'chemotherapy', 'chemokine']

        self.parameters = params
        self.attributes = atts

        blood_vessels = [[0, 0]]
        fast_bacteria = [[1, 0], [1, 1]]
        slow_bacteria = [[2, 0], [2, 1]]
        macrophages = [[3, 0], [3, 1], [3, 2], [3, 3]]
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [12, 12], atts, params, blood_vessels, fast_bacteria,
                                                        slow_bacteria, macrophages)

        self.topology.automata[0].grid[1, 1]['contents'].resting = True
        self.topology.automata[0].grid[2, 1]['contents'].resting = True

        self.topology.automata[0].grid[3, 1]['contents'].state = 'active'
        self.topology.automata[0].grid[3, 2]['contents'].state = 'infected'
        self.topology.automata[0].grid[3, 3]['contents'].state = 'chronically_infected'

        tcell = TB_Model.TCell([4,0])
        self.topology.automata[0].t_cells.append(tcell)
        self.topology.automata[0].grid[4,0]['contents'] = tcell

        self.topology.automata[0].caseum.append([5,0])
        self.topology.automata[0].grid[5,0]['contents'] = 'caseum'

        self.topology.automata[0].max_oxygen_global = 1.5

    def sort_out_halos(self):
        dz = []
        for i in self.topology.automata:
            dz.append(i.get_danger_zone())
        halos = self.topology.create_halos(dz)
        for i in range(4):
            self.topology.automata[i].set_halo(halos[i])

    def test_output(self):

        self.topology.automata[0].record_state()

        # Check the files exist
        os.path.isfile('0_contents.txt')
        os.path.isfile('0_oxygen.txt')
        os.path.isfile('0_chemotherapy.txt')
        os.path.isfile('0_chemokine.txt')

        # Check contents file
        content_lines = [line.rstrip('\n') for line in open('0_contents.txt')]
        for i in range(len(content_lines)):

            if i == 6:
                self.assertEqual(content_lines[i], '1.0')
            elif i == 7:
                self.assertEqual(content_lines[i], '1.25')
            elif i == 12:
                self.assertEqual(content_lines[i], '2.0')
            elif i == 13:
                self.assertEqual(content_lines[i], '2.25')
            elif i == 18:
                self.assertEqual(content_lines[i], '4.0')
            elif i == 19:
                self.assertEqual(content_lines[i], '5.0')
            elif i == 20:
                self.assertEqual(content_lines[i], '6.0')
            elif i == 21:
                self.assertEqual(content_lines[i], '7.0')
            elif i == 24:
                self.assertEqual(content_lines[i], '3.0')
            elif i == 30:
                self.assertEqual(content_lines[i], '100.0')
            else:
                self.assertEqual(content_lines[i], '0.0')

        # Check oxygen file
        oxygen_lines = [line.rstrip('\n') for line in open('0_oxygen.txt')]
        for i in range(len(oxygen_lines)):
            if i == 0:
                self.assertEqual(oxygen_lines[i], '100.0')
            else:
                self.assertEqual(oxygen_lines[i], '0.0')

        # CHEMOTHERAPY
        chemotherapy_lines = [line.rstrip('\n') for line in open('0_chemotherapy.txt')]
        for i in range(len(chemotherapy_lines)):
            self.assertEqual(chemotherapy_lines[i], '0.0')

        # CHEMOKINE
        chemokine_lines = [line.rstrip('\n') for line in open('0_chemokine.txt')]
        for i in range(len(chemokine_lines)):
            self.assertEqual(chemokine_lines[i], '0.0')

        self.sort_out_halos()
        self.topology.automata[0].update()
        self.topology.automata[0].process_events([])

    def test_parameter(self):

        self.sort_out_halos()
        self.topology.automata[0].parameters['interval_to_record_results'] = 2

        self.topology.automata[0].update()
        self.topology.automata[0].process_events([])

        self.assertFalse(os.path.isfile('0_contents.txt'))
        self.assertFalse(os.path.isfile('0_oxygen.txt'))
        self.assertFalse(os.path.isfile('0_chemotherapy.txt'))
        self.assertFalse(os.path.isfile('0_chemokine.txt'))

        self.topology.automata[0].update()
        self.topology.automata[0].process_events([])

        self.assertTrue(os.path.isfile('0_contents.txt'))
        self.assertTrue(os.path.isfile('0_oxygen.txt'))
        self.assertTrue(os.path.isfile('0_chemotherapy.txt'))
        self.assertTrue(os.path.isfile('0_chemokine.txt'))

if __name__ == '__main__':
    unittest.main()
