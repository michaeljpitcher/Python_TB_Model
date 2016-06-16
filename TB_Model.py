import numpy as np
import math

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