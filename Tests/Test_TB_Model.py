import unittest
import TB_Model


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

    def test_address_to_integer(self):
        x = 0
        y = 0
        for i in range(25):
            self.assertEqual(self.tile.address_to_location([x, y]), i)
            y += 1
            if y == 5:
                y = 0
                x += 1

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
        self.tile.set_attribute_work_grid([4, 4], 'a', 99)
        self.tile.set_attribute_work_grid([4, 3], 'b', 99)
        self.tile.set_attribute_work_grid([3, 4], 'c', 99)

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
        # TODO - expect exception
        # self.tile.set_attribute_grid([2, 2], 'not set', 99.0)

    def test_set_att_work(self):
        self.tile.set_attribute_work_grid([2, 2], 'a', 99.0)
        self.assertEqual(self.tile.work_grid[2, 2]['a'], 99.0)
        # TODO - expect exception
        # self.tile.set_attribute_grid([2, 2], 'not set', 99.0)


class NeighbourhoodTestCase(unittest.TestCase):

    def setUp(self):
        self.neighbourhood = TB_Model.Neighbourhood(2, 3)

    def test_neighbour_table(self):
        table = self.neighbourhood.neighbour_table[1]
        self.assertItemsEqual(table, [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)])

    def test_calculate_neighbour_location(self):
        table = self.neighbourhood.neighbour_table[1]
        locations = self.neighbourhood.calculate_neighbours_locations([2, 2], table)
        self.assertItemsEqual(locations, [[1, 1], [1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2], [3, 3]])

    def test_neighbours_moore(self):
        nm_1 = self.neighbourhood.neighbours_moore([2, 2], 1, True)
        self.assertItemsEqual(nm_1, [[1, 1], [1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2], [3, 3]])
        nm_2_inc = self.neighbourhood.neighbours_moore([2, 2], 2, True)
        self.assertItemsEqual(nm_2_inc, [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4],
                                         [1, 0], [1, 1], [1, 2], [1, 3], [1, 4],
                                         [2, 0], [2, 1], [2, 3], [2, 4],
                                         [3, 0], [3, 1], [3, 2], [3, 3], [3, 4],
                                         [4, 0], [4, 1], [4, 2], [4, 3], [4, 4]])
        nm_2_exc = self.neighbourhood.neighbours_moore([2, 2], 2, False)
        self.assertItemsEqual(nm_2_exc, [[0, 0], [0, 1], [0, 2], [0, 3], [0, 4],
                                         [1, 0], [1, 4],
                                         [2, 0], [2, 4],
                                         [3, 0], [3, 4],
                                         [4, 0], [4, 1], [4, 2], [4, 3], [4, 4]])

    def test_neighbours_von_neumann(self):
        nvn_1 = self.neighbourhood.neighbours_von_neumann([2, 2], 1, True)
        self.assertItemsEqual(nvn_1, [[1, 2], [2, 1], [2, 3], [3, 2]])
        nvn_2_inc = self.neighbourhood.neighbours_von_neumann([2, 2], 2, True)
        self.assertItemsEqual(nvn_2_inc,
                              [[0, 2], [1, 1], [1, 2], [1, 3], [2, 0], [2, 1], [2, 3], [2, 4], [3, 1], [3, 2], [3, 3],
                               [4, 2]])
        nvn_2_exc = self.neighbourhood.neighbours_von_neumann([2, 2], 2, False)
        self.assertItemsEqual(nvn_2_exc, [[0, 2], [1, 1], [1, 3], [2, 0], [2, 4], [3, 1], [3, 3], [4, 2]])


class AutomatonTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        atts = ['a','b','c']
        self.automaton = TB_Model.Automaton([5,5], 1, atts, params)


class TopologyTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        atts = ['a', 'b', 'c', 'blood_vessel']
        self.topology = TB_Model.Topology([2,2],[10,10],atts,params,[[[3,3]],[],[],[]])

    def test_init(self):
        self.assertEqual(len(self.topology.automata),4)
        self.assertEqual(self.topology.automata[0].shape[0], 5)
        self.assertEqual(self.topology.automata[0].shape[1], 5)
        self.assertEqual(self.topology.automata[1].shape[0], 5)
        self.assertEqual(self.topology.automata[1].shape[1], 5)
        self.assertEqual(self.topology.automata[2].shape[0], 5)
        self.assertEqual(self.topology.automata[2].shape[1], 5)
        self.assertEqual(self.topology.automata[3].shape[0], 5)
        self.assertEqual(self.topology.automata[3].shape[1], 5)

    def test_get_external_addresses_required(self):
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
        gear = self.topology.get_external_addresses_required(self.topology.automata[0])
        self.assertItemsEqual(addresses, gear)


class TwoDimensionalTopologyTestCase(unittest.TestCase):

    def setUp(self):
        params = dict()
        params['max_depth'] = 3
        atts = ['a', 'b', 'c', 'blood_vessel']
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, [[3,3]])

    def test_init(self):
        self.assertSequenceEqual(self.topology.origins,[[0,0], [0,5], [5,0], [5,5]])

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
        self.attributes = ['a', 'blood_vessel']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters,[[3,3]])

        self.assertSequenceEqual(self.topology.local_to_local(0, [0, 4], 1), [0, 1])
        self.assertSequenceEqual(self.topology.local_to_local(0, [3, 2], 2), [0, 2])
        self.assertSequenceEqual(self.topology.local_to_local(0, [3, 3], 3), [0, 0])

        self.assertSequenceEqual(self.topology.local_to_local(0, [0, 2], 1), [0, -1])
        self.assertSequenceEqual(self.topology.local_to_local(1, [0, 0], 2), [-3, 3])

        self.assertSequenceEqual(self.topology.local_to_local(3, [-1, -1], 0), [2, 2])
        self.assertSequenceEqual(self.topology.local_to_local(3, [-1, 0], 1), [2, 0])
        self.assertSequenceEqual(self.topology.local_to_local(3, [0, -2], 2), [0, 1])

    def test_create_halos(self):
        self.attributes = ['a', 'blood_vessel']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters, [[3,3]])

        tile_id = 0
        x = 0
        y = 0
        for i in range(36):
            self.topology.automata[tile_id].set_attribute_work_grid([x, y], 'a', i)
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
        self.attributes = ['a', 'blood_vessel']
        self.parameters = dict()
        self.parameters['max_depth'] = 1
        self.topology = TB_Model.TwoDimensionalTopology([2, 2], [6, 6], self.attributes, self.parameters, [[3,3]])

        tile_id = 0
        x = 0
        y = 0
        for i in range(36):
            self.topology.automata[tile_id].set_attribute_work_grid([x, y], 'a', i)
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


if __name__ == '__main__':
    unittest.main()
