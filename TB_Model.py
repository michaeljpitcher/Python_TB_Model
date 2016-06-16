import numpy as np
import math
import itertools


class Topology:

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local = None):

        self.number_of_tiles = reduce(lambda x, y: x * y, tile_arrangement)
        self.total_shape = np.array(total_shape)
        self.tile_shape = self.total_shape / tile_arrangement
        self.tile_arrangement = tile_arrangement

        self.automata = []

        for i in range(self.number_of_tiles):
            # TODO
            automaton = Automaton(self.tile_shape, i, attributes, parameters, blood_vessel_local[i], [], [], [])
            self.automata.append(automaton)

        self.external_addresses_required = self.get_external_addresses_required(self.automata[0])

        for automaton in self.automata:
            automaton.halo_addresses = self.external_addresses_required

    def get_external_addresses_required(self, automaton):
        """
            The addresses of cells within the overall neighbourhood which don't belong to this automaton
            :return:
            """
        external_addresses = []
        # Loop through every cell
        for i in range(0, automaton.grid.size):
            # Pull the addresses of cells in neighbourhood of this cell
            address = automaton.location_to_address(i)
            neighbours = automaton.neighbours_moore(address, automaton.parameters['max_depth'])
            # For every neighbour cell
            for neighbour in neighbours:
                # Check neighbour not on grid and hasn't been processes already
                if (not automaton.address_is_on_grid(neighbour)) and neighbour not in external_addresses:
                    external_addresses.append(neighbour)

        return external_addresses


class TwoDimensionalTopology(Topology):
    """
    2d Grid Topology (amend normalise address if Toroid needed)
    """

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_global = None):
        assert len(total_shape) == 2
        self.number_of_tiles = reduce(lambda x, y: x * y, tile_arrangement)
        self.total_shape = np.array(total_shape)
        self.tile_shape = self.total_shape / tile_arrangement
        self.tile_arrangement = tile_arrangement
        blood_vessel_local = self.get_blood_vessel_local(blood_vessel_global)
        Topology.__init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local)

        # Create a list detailing where each tile's origin (local 0,0) lies in relation to the global grid
        self.origins = []
        x = 0
        y = 0
        for z in range(self.number_of_tiles):
            self.origins.append([x, y])
            y += self.tile_shape[1]
            # If we've reached the width of the grid, reset y to 0 and increment x (start a new row)
            if y % self.total_shape[1] == 0:
                y = 0
                x += self.tile_shape[0]

        # Get a list of every global address needed
        self.global_addresses_required = []
        for automaton in self.automata:
            for b in automaton.halo_addresses:
                address = self.local_to_global(automaton.tile_id, b)
                if address is not None:
                    self.global_addresses_required.append(address)

        # Use required global addresses to create danger zones
        self.danger_zone_addresses = dict()
        for tile_id in range(self.number_of_tiles):
            self.danger_zone_addresses[tile_id] = []
        for global_address in self.global_addresses_required:
            tile_id, local_address = self.global_to_local(global_address)
            if local_address not in self.danger_zone_addresses[tile_id]:
                self.danger_zone_addresses[tile_id].append(local_address)

        # Set the danger zone addresses
        for tile_id in self.danger_zone_addresses.keys():
            self.automata[tile_id].set_addresses_for_danger_zone(self.danger_zone_addresses[tile_id])

    def normalise_address(self, address):
        """
        Normalise the address
        :param address:
        :return: None if outside the global boundary, else the address
        """
        # Normalise the address - returns None if outside the global boundary
        x, y = address
        if x < 0 or x >= self.total_shape[0] or y < 0 or y >= self.total_shape[1]:
            return None
        return [x, y]

    def global_to_local(self, global_address):
        if global_address is None:
            return None

        x, y = global_address
        tile_rows, tile_cols = self.tile_shape

        # Add the tile ID - x modulo num of rows in a tile * number of tiles in width of the grid
        output = [divmod(x, tile_rows)[0] * self.tile_arrangement[1] + divmod(y, tile_cols)[0]]
        # Add the Local coordinates
        output.append([divmod(x, tile_rows)[1], divmod(y, tile_cols)[1]])

        return output

    def local_to_global(self, tile_id, local_address):
        x, y = local_address
        origin = self.origins[tile_id]
        # Normalise the address before it's returned
        return self.normalise_address([origin[0] + x, origin[1] + y])

    def local_to_local(self, original_tile_id, local_address, new_tile_id):
        global_address = self.local_to_global(original_tile_id, local_address)
        origin_x, origin_y = self.origins[new_tile_id]
        return [global_address[0] - origin_x, global_address[1] - origin_y]

    def create_halos(self, danger_zone_values):
        """
        Get danger zones, turn to halos
        :return:
        """

        global_values_required = []

        for g in self.global_addresses_required:
            tile_id, local_address = self.global_to_local(g)
            index = self.danger_zone_addresses[tile_id].index(local_address)
            cell = danger_zone_values[tile_id][index]
            global_values_required.append(cell)

        halos = []

        for tile_id in range(self.number_of_tiles):
            halo = []
            for address_required in self.external_addresses_required:
                global_address = self.local_to_global(tile_id, address_required)
                if global_address is None:
                    halo.append(None)
                else:
                    index = self.global_addresses_required.index(global_address)
                    value = global_values_required[index]
                    halo.append(value)
            halos.append(halo)

        return halos

    def get_blood_vessel_local(self, blood_vessel_global):

        blood_vessel_local = []
        for i in range(self.number_of_tiles):
            blood_vessel_local.append([])

        for global_address in blood_vessel_global:
            tile_id, local_address = self.global_to_local(global_address)
            blood_vessel_local[tile_id].append(local_address)

        return blood_vessel_local


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

    # def initialise(self, attributes, values):
    #
    #     # Set all the initial values
    #     for index in range(len(values)):
    #         value_list = values[index]
    #         assert len(value_list) == self.size
    #         attribute = attributes[index]
    #         for cell_index in range(len(value_list)):
    #             value = value_list[cell_index]
    #             address = self.location_to_address(cell_index)
    #             self.set_attribute_grid(address, attribute, value)
    #
    #     # Clone the work grid
    #     self.create_work_grid()


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

    def __init__(self, shape, tile_id, attributes, parameters, blood_vessels, fast_bacteria=None, slow_bacteria=None,
                 macrophages=None):
        Tile.__init__(self, shape, attributes)
        Neighbourhood.__init__(self, len(shape), parameters['max_depth'])
        self.tile_id = tile_id
        self.parameters = parameters

        # INITIAL
        self.initialise_blood_vessels(blood_vessels)
        # TODO - agents initial

        # COPY GRID TO WORK GRID
        self.create_work_grid()

    def initialise_blood_vessels(self, addresses):

        for address in addresses:
            self.set_attribute_grid(address,'blood_vessel',1.5)




