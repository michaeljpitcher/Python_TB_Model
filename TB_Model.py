import numpy as np
import math
import itertools

class Tile:

    def __init__(self, shape, attributes):

        self.shape = shape
        self.attributes = attributes
        self.size = reduce(lambda x, y: x * y, shape)

        self.create_grid(attributes)

    def create_grid(self, attributes):

        cells = []
        for i in range(self.size):
            cell = dict()
            for att in attributes:
                cell[att] = 0.0
            cells.append(cell)

        self.grid = np.array(cells).reshape(self.shape)

    def create_work_grid(self):

        cells = []
        for i in range(self.size):
            address = self.location_to_address(i)
            cell = self.grid[tuple(address)].copy()
            cells.append(cell)

        self.work_grid = np.array(cells).reshape(self.shape)

    def swap_grids(self):
        """
        Swap the active grid for the working grid
        :return:
        """
        self.grid, self.work_grid = self.work_grid, self.grid

    def location_to_address(self, integer):
        """
        Turn an integer value into an address list of co-ordinates
        :param integer:
        :return:
        """
        return list(np.unravel_index(integer, self.grid.shape))

    def address_to_location(self, address):
        """
        Convert an address set of coordinates into an integer location
        :param address:
        :return:
        """
        result = 0
        acc = 1
        for pi, si in zip(reversed(address), reversed(self.grid.shape)):
            result += pi * acc
            acc *= si
        return result

    def address_is_on_grid(self, address):
        """
        Check if address is within the boundaries
        :param address:
        :return:
        """
        for i in range(len(address)):
            if address[i] < 0 or address[i] >= self.shape[0]:
                return False

        return True

    def set_addresses_for_danger_zone(self, addresses):
        self.danger_zone_addresses = addresses

    def get_danger_zone(self):
        danger_zone = []
        for address in self.danger_zone_addresses:
            # Get cell from the work grid
            danger_zone.append(self.work_grid[tuple(address)])

        return danger_zone

    def set_addresses_for_halo(self, addresses):
        self.halo_addresses = addresses

    def set_halo(self, cells):
        self.halo_cells = []
        for i in range(len(cells)):
            self.halo_cells.append(cells[i])

    def get(self, address):
        """
        Get a cell from the active grid
        :param address:
        :return:
        """
        if self.address_is_on_grid(address):
            address = tuple(address)
            return self.grid[address]
        elif address in self.halo_addresses:
            index = self.halo_addresses.index(address)
            return self.halo_cells[index]
        else:
            raise Exception("Failure at get method")

    def get_attribute(self, address, attribute):
        """
        Get an attribute from a cell on the active grid
        :param address:
        :param attribute:
        :return:
        """
        if self.address_is_on_grid(address):
            address = tuple(address)
            return self.grid[address][attribute]
        elif address in self.halo_addresses:
            index = self.halo_addresses.index(address)
            if self.halo_cells[index] is None:
                return None
            else:
                return self.halo_cells[index][attribute]
        else:
            raise Exception("Failure at get method")

    def set_attribute_grid(self, address, attribute, value):
        """
        Set an attribute of a cell on the active grid (for initialisation)
        :param address:
        :param attribute:
        :param value:
        :return:
        """
        address = tuple(address)
        if attribute in self.grid[address].keys():
            self.grid[address][attribute] = value
        else:  # Specified attribute hasn't been set as a possibility
            raise Exception('Attribute {0} does not exist'.format(attribute))

    def set_attribute_work_grid(self, address, attribute, value):
        """
        Set an attribute of a cell on the working grid
        :param address:
        :param attribute:
        :param value:
        :return:
        """
        address = tuple(address)
        if attribute in self.work_grid[address].keys():
            self.work_grid[address][attribute] = value
        else:  # Specified attribute hasn't been set as a possibility
            raise Exception('Attribute {0} does not exist'.format(attribute))

    def initialise(self, attributes, values):

        # Set all the initial values
        for index in range(len(values)):
            value_list = values[index]
            assert len(value_list) == self.size
            attribute = attributes[index]
            for cell_index in range(len(value_list)):
                value = value_list[cell_index]
                address = self.location_to_address(cell_index)
                self.set_attribute_grid(address, attribute, value)

        # Clone the work grid
        self.create_work_grid()

class Neighbourhood:

    def __init__(self, dimensions, max_depth):

        self.dimensions = dimensions
        self.neighbour_table = self.construct_neighbour_table(max_depth)

    def construct_neighbour_table(self, max_depth):
        """
        Create a list of relative neighbour addresses
        :param depth:
        :return:
        """

        table = dict()

        for d in range(max_depth):

            depth = d+1
            range_ = range(-depth, depth + 1)
            row = list(itertools.product(range_, repeat=self.dimensions))
            row.remove((0,) * self.dimensions)
            table[depth] = row

        return table

    def calculate_neighbours_locations(self, address, table):
        """
        Applies the given neighbour table to the given address to get addresses of neighbours
        :param address:
        :param table:
        :return:
        """
        output = []
        for i in range(len(table)):
            # Build new neighbour
            neighbour = []
            for j in range(self.dimensions):
                neighbour.append(address[j] + table[i][j])
            output.append(neighbour)
        return output

    def neighbours_moore(self, address, depth=1, inclusive=True):
        """
        Address of neighbours in the Moore neighbourhood of given depth
        :param address:
        :param depth:
        :param inclusive: Include lower depths?
        :return:
        """
        table = self.neighbour_table[depth]
        # If not inclusive, remove any neighbours in smaller depths
        # TODO - is this a slow way of doing things?
        if not inclusive:
            reduced_table = []
            for neighbour in table:
                for x in neighbour:
                    if int(math.fabs(x)) == depth:
                        reduced_table.append(neighbour)
                        break
            table = reduced_table
        return self.calculate_neighbours_locations(address, table)

    def neighbours_von_neumann(self, address, depth=1, inclusive=True):
        """
        Gives the neighbours in the Von Neumann neighbourhood of address of given depth
        :param address:
        :param depth:
        :param inclusive: Includes lower depths?
        :return:
        """

        # Build truth table of all possible values based on depth of the neighbourhood and number of dimensions
        table = self.neighbour_table[depth]
        # Reduce the values based on their Manhattan distance
        # TODO - is this a slow way of doing things?
        reduced_table = []
        for k in table:
            # Check the Manhattan distance is up to the depth (inclusive) or equal to depth (exclusive)
            if (inclusive and int(sum([math.fabs(x) for x in k])) <= depth) or \
                    ((not inclusive) and int(sum([math.fabs(x) for x in k])) == depth):
                reduced_table.append(k)

        return self.calculate_neighbours_locations(address, reduced_table)


class Automaton(Tile, Neighbourhood):

    def __init__(self, shape, tile_id, attributes, parameters):
        Tile.__init__(self, shape, attributes)
        Neighbourhood.__init__(self, len(shape), parameters['max_depth'])
        self.tile_id = tile_id
        self.parameters = parameters





