import unittest
import TB_Model


class TileTestCase(unittest.TestCase):

    def setUp(self):
        self.tile = TB_Model.Tile([5, 5], ['a', 'b', 'c'])

    def test_shape(self):
        self.assertEqual(self.tile.grid.shape, (5, 5))

    def test_create_work_grid(self):
        self.tile.create_work_grid()
        
        self.tile.grid[0,0]['a'] = 1.1
        print self.tile.work_grid[0,0]

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
        self.tile.set_attribute_work_grid([0, 0], 'a', 1)
        self.tile.swap_grids()
        self.assertEqual(self.tile.grid[0, 0]['a'], 1)
        self.assertEqual(self.tile.work_grid[0, 0]['a'], 0.0)

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

    def test_set_DZ(self):

        self.tile.set_addresses_for_danger_zone([[4, 4], [4, 3], [3, 4]])
        danger_zone = self.tile.get_danger_zone()

        for i in danger_zone:
            i['c'] = 77

        self.tile.set_danger_zone(danger_zone)

        self.assertEqual(self.tile.work_grid[4, 4]['c'], 77)
        self.assertEqual(self.tile.work_grid[4, 3]['c'], 77)
        self.assertEqual(self.tile.work_grid[3, 4]['c'], 77)

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
        self.assertItemsEqual(self.tile.get([0, 0]).keys(), ['a', 'b', 'c'])
        self.assertItemsEqual(self.tile.get([0, 0]).values(), [0.0, 0.0, 0.0])

    def test_get_attribute_on_grid(self):
        self.assertEqual(self.tile.get_attribute([0, 0], 'a'), 0.0)

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

    def test_set(self):
        cell = dict()
        cell['d'] = 1
        cell['e'] = 2
        self.tile.set_grid([2, 2], cell)
        self.assertItemsEqual(self.tile.grid[2, 2].keys(), ['d', 'e'])
        self.assertEqual(self.tile.grid[2, 2]['d'], 1)
        self.assertEqual(self.tile.grid[2, 2]['e'], 2)

    def test_set_work(self):
        cell = dict()
        cell['d'] = 1
        cell['e'] = 2
        self.tile.set_work_grid([2, 2], cell)
        self.assertItemsEqual(self.tile.work_grid[2, 2].keys(), ['d', 'e'])
        self.assertEqual(self.tile.work_grid[2, 2]['d'], 1)
        self.assertEqual(self.tile.work_grid[2, 2]['e'], 2)

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



if __name__ == '__main__':
    unittest.main()
