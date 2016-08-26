import numpy as np
import math
import itertools
from collections import Counter
import os

'''
Tuberculosis Automaton Model
Michael J. Pitcher
School of Computer Science, University of St. Andrews

Creates a parallelised cellular automaton / agent hybrid model, where the grid is split into multiple smaller grids
to improve performance. A series of automata run updates on the tiles and agents which are on the grids, and these
create a series of potential events. Once conflicting events are resolved (handled elsewhere) the acceptable events
are passed back to automata which update the tiles, and the process repeats as needed.

Some terminology:

Cell - the individual blocks of the tile (here represented as a dictionary of attributes)
Tile - the smaller grids which constitute the larger overall grid
Size - Number of cells in a tile
Shape - Arrangement of cells in a tile
Address - an (x,y,z,etc.) collection of co-ordinates for a cell
Halo - the cells required by a tile for updating which are not a part of the tile (belong to other tiles)
Danger Zone - the cells in a tile which will be required by other tiles (i.e. are part of other tiles' halos)
'''


class Topology:
    """
    Class that handles data transfer between a set of distinct automata. Determines what cells each tile needs from
    other tiles and constructs the halos that are passed between
    """

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local,
                 fast_bacteria_local, slow_bacteria_local, macrophages_local):

        self.number_of_tiles = reduce(lambda x, y: x * y, tile_arrangement)
        self.total_shape = np.array(total_shape)
        self.tile_shape = self.total_shape / tile_arrangement
        self.tile_arrangement = tile_arrangement

        # List of automata - each tile acts as it's own automata and runs updates independently of other tiles
        self.automata = []

        # Create automata
        for i in range(self.number_of_tiles):
            automaton = Automaton(self.tile_shape, i, attributes, parameters, blood_vessel_local[i],
                                  fast_bacteria_local[i], slow_bacteria_local[i], macrophages_local[i])
            self.automata.append(automaton)

        # Get a list of addresses (relative to tile) that are required by tile but outside of it
        # Will be the same for all tiles (as they're same shape) so just pull from one
        self.external_addresses_required = self.get_external_addresses_required(parameters['max_depth'])
        # Halo of depth 1 - needed for calculating diffusion rates
        # TODO - COMP - is there a better way of doing this?
        depth1_addresses = self.get_external_addresses_required(1)

        # Set the halo address and halo of depth 1 address
        for automaton in self.automata:
            automaton.configure_halo_addresses(self.external_addresses_required, depth1_addresses)
            automaton.configure_neighbourhood_for_halo(self.external_addresses_required)

    def get_external_addresses_required(self, depth):
        """
        The addresses of cells within the overall neighbourhood which don't belong to the given automaton
        :return:
        """
        external_addresses = []
        automaton = self.automata[0]
        # Loop through every cell
        for address in automaton.list_grid_addresses:
            # Pull the addresses of cells in neighbourhood of this cell
            neighbours = []
            for j in range(1, depth+1):
                neighbours += automaton.moore_neighbours[address][j]
            # For every neighbour cell
            for neighbour in neighbours:
                # Check neighbour not on grid and hasn't been processes already
                if (not automaton.address_is_on_grid(neighbour)) and neighbour not in external_addresses:
                    external_addresses.append(neighbour)

        return external_addresses


class TwoDimensionalTopology(Topology):
    """
    2d Grid Topology
    """

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_global,
                 fast_bacteria_global, slow_bacteria_global, macrophages_global):

        for i in range(len(total_shape)):
            if total_shape[i] % tile_arrangement[i] != 0:
                raise Exception("Dimension {0} does not divide evenly by tile arrangement {1}".format(total_shape[i],
                                                                                                      tile_arrangement[
                                                                                                          i]))
        # Check it's two dimensions
        assert len(total_shape) == 2
        self.number_of_tiles = reduce(lambda a, q: a * q, tile_arrangement)
        self.total_shape = np.array(total_shape)
        # Calculate each tile shape by dividing total shape by tile arrangement
        self.tile_shape = self.total_shape / tile_arrangement
        self.tile_arrangement = tile_arrangement

        # Initialise (turn the global addresses into lists of local ones)
        blood_vessel_local = self.get_local_addresses(blood_vessel_global)
        fast_bacteria_local = self.get_local_addresses(fast_bacteria_global)
        slow_bacteria_local = self.get_local_addresses(slow_bacteria_global)
        macrophages_local = self.get_local_addresses(macrophages_global)

        Topology.__init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local,
                          fast_bacteria_local, slow_bacteria_local, macrophages_local)

        # Create a list detailing where each tile's origin (local 0,0) lies in relation to the global grid
        self.origins = []
        x = 0
        y = 0
        for z in range(self.number_of_tiles):
            self.origins.append((x, y))
            y += self.tile_shape[1]
            # If we've reached the width of the grid, reset y to 0 and increment x (start a new row)
            if y % self.total_shape[1] == 0:
                y = 0
                x += self.tile_shape[0]

        # Get a list of every global address needed - every address that's on one tile but needed by another tile
        self.global_addresses_required = []
        for automaton in self.automata:
            # for b in automaton.halo_addresses:
            for b in automaton.list_halo_addresses:
                address = self.local_to_global(automaton.tile_id, b)
                if address is not None:
                    self.global_addresses_required.append(address)

        # Use required global addresses to create danger zones on tiles
        # Danger zones are addresses on a tile that are required by other tiles
        self.danger_zone_addresses = dict()
        # Set up empty lists
        for tile_id in range(self.number_of_tiles):
            self.danger_zone_addresses[tile_id] = []
        # Loop through each of the global addresses required, convert them into local coordinates
        for global_address in self.global_addresses_required:
            tile_id, local_address = self.global_to_local(global_address)
            if local_address not in self.danger_zone_addresses[tile_id]:
                self.danger_zone_addresses[tile_id].append(local_address)

        # Set the danger zone addresses
        for tile_id in self.danger_zone_addresses.keys():
            self.automata[tile_id].set_addresses_for_danger_zone(self.danger_zone_addresses[tile_id])

    def global_to_local(self, global_address):
        """
        Turn a global address into a tile ID and a local address
        :param global_address:
        :return:
        """
        if global_address is None:
            return None, None

        x, y = global_address
        tile_rows, tile_cols = self.tile_shape

        # Add the tile ID - x modulo num of rows in a tile * number of tiles in width of the grid
        output = [divmod(x, tile_rows)[0] * self.tile_arrangement[1] + divmod(y, tile_cols)[0]]
        # Add the Local coordinates
        x = divmod(x, tile_rows)[1]
        y = divmod(y, tile_cols)[1]
        output.append((x, y))

        return output

    def local_to_global(self, tile_id, local_address):
        """
        Turn a tile ID and a local address into a global address
        :param tile_id:
        :param local_address:
        :return:
        """
        origin = self.origins[tile_id]

        new_address = (origin[0] + local_address[0], origin[1] + local_address[1])
        # Normalise the address before it's returned
        if new_address[0] < 0 or new_address[0] >= self.total_shape[0] or \
                new_address[1] < 0 or new_address[1] >= self.total_shape[1]:
            return None
        return new_address

    def local_to_local(self, original_tile_id, local_address, new_tile_id):
        """
        Turn a tile ID and a local address into a local address for the given tile
        Basically does local -> global -> local
        :param original_tile_id:
        :param local_address:
        :param new_tile_id:
        :return:
        """
        global_address = self.local_to_global(original_tile_id, local_address)
        origin_x, origin_y = self.origins[new_tile_id]
        return global_address[0] - origin_x, global_address[1] - origin_y

    def create_halos(self, danger_zone_values):
        """
        Turn the danger zone cell values into halo values. Uses the pre-determined addresses.
        :param danger_zone_values: The values of cells in danger zone.
        :return:
        """
        global_cells_required = dict()

        # From the lists of danger zone values, turn each into an entry in the global cells required dictionary
        # Check DZ cells for each tile
        for tile_id in range(self.number_of_tiles):
            # Get the addresses for this tile
            dz_addresses = self.danger_zone_addresses[tile_id]
            dz_cells = danger_zone_values[tile_id]
            for index in range(len(dz_cells)):
                # Convert to a global coordinate and store
                global_address = self.local_to_global(tile_id, dz_addresses[index])
                global_cells_required[global_address] = dz_cells[index]

        # Now have a list of global addresses and their cell values. Turn these into the halos
        halos = []
        for tile_id in range(self.number_of_tiles):
            halo = []
            # Loop through each required halo address
            for address_required in self.external_addresses_required:
                # Find it's global entry and store
                global_address = self.local_to_global(tile_id, address_required)
                # Store none if off the grid
                if global_address is None:
                    halo.append(None)
                else:
                    value = global_cells_required[global_address]
                    halo.append(value)
            halos.append(halo)
        return halos

    def get_local_addresses(self, global_addresses):
        """
        Turn a set of global address into local addresses
        :param global_addresses:
        :return:
        """

        local_addresses = []
        for i in range(self.number_of_tiles):
            local_addresses.append([])

        for global_address in global_addresses:
            tile_id, local_address = self.global_to_local(global_address)
            local_addresses[tile_id].append(local_address)

        return local_addresses


class Tile:
    """
    The data storage of the automata. Contains the grid where values and events are calculated from, the work grid where
    new values are entered into
    """

    def __init__(self, shape, attributes):
        self.shape = shape
        self.attributes = attributes
        self.size = reduce(lambda x, y: x * y, shape)

        # List of addresses - stored to allow easy iteration through every cell
        self.list_grid_addresses = []
        self.list_halo_addresses = []

        # Grids are dictionaries, keys are coordinates as tuples e.g. (0,0) and values are dictionaries of attributes
        self.grid = dict()
        self.work_grid = dict()
        # Populate address list and grid
        for i in range(self.size):
            # Turn integer value into address tuple
            address = np.unravel_index(i, self.shape)
            # Cell value is a dictionary, where each key is an attribute
            self.grid[address] = dict()
            for a in attributes:
                self.grid[address][a] = 0.0
            self.list_grid_addresses.append(address)
        self.danger_zone_addresses = []
        self.halo_depth1 = []

    def create_work_grid(self):
        """
        Create a working grid placeholder by cloning the current grid
        :return:
        """
        work_grid = dict()
        for address in self.list_grid_addresses:
            # Need to copy() else just points to original grid
            cell = self.grid[address].copy()
            work_grid[address] = cell
        self.work_grid = work_grid

    def swap_grids(self):
        """
        Swap the active grid for the working grid
        :return:
        """
        self.grid, self.work_grid = self.work_grid, self.grid

    def address_is_on_grid(self, address):
        """
        Check if address is within the boundaries of the tile
        :param address:
        :return:
        """
        # Compare each coordinate against the shape, if any one is outside return False
        for i in range(len(address)):
            if address[i] < 0 or address[i] >= self.shape[i]:
                return False
        return True

    def set_addresses_for_danger_zone(self, addresses):
        self.danger_zone_addresses = addresses

    def get_danger_zone(self):
        """
        Get the cell values of cells in the danger zone
        :return:
        """
        # Pull each address in stored danger zone addresses and get its value from grid
        return [self.grid[address] for address in self.danger_zone_addresses]

    def configure_halo_addresses(self, external_addresses_required, depth1_addresses):
        """
        Given halo addresses, add them to the grid structure for easy reference
        :param external_addresses_required:
        :param depth1_addresses:
        :return:
        """
        # By adding halo to the grid, means don't need to check for where a cell is (i.e. if its on tile or in halo)
        # For each halo address, add an empty cell to grid and add address to list
        for halo_address in external_addresses_required:
            self.grid[halo_address] = dict()
            self.list_halo_addresses.append(halo_address)
        # Store the depth 1 halo
        self.halo_depth1 = depth1_addresses

    def set_halo(self, halo):
        """
        Given a set of halo cells, assign the cells to the grid in correct locations
        :param halo:
        :return:
        """
        # Use the already stored list_halo_addresses as guide for where to store each cell value
        for index in range(len(halo)):
            address = self.list_halo_addresses[index]
            self.grid[address] = halo[index]

    def set_attribute_grid(self, address, attribute, value):
        """
        Set an attribute of a cell on the active grid (for initialisation and pre-processing)
        :param address:
        :param attribute:
        :param value:
        :return:
        """
        # Check that the attribute exists, if so, then set the value to given cell
        if attribute in self.attributes:
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
        self.work_grid[address][attribute] = value


class Neighbourhood:
    """
    Data structure to calculate addresses. Upon initialisation, loops through every cell and calculates all neighbours
    that may be needed and stores them in a dictionary for future lookup.
    """

    def __init__(self, dimensions, max_depth, list_addresses):
        self.max_depth = max_depth
        self.dimensions = dimensions

        # Dictionaries of neighbour addresses. Keys are addresses of cells, values are dictionaries - these have keys
        # 1 to max_depth and values are lists of neighbours
        # e.g. Moore neighbours of address (1,1) to a depth of 2 would be self.moore_neighbours[(1,1)][2]
        # THESE ARE NON-INCLUSIVE (to save data size and make searching for space at succesive levels easier), so a full
        # moore neighbourhood of depth 2 needs to look at depth 1 and then depth 2
        self.moore_neighbours = dict()
        self.von_neumann_neighbours = dict()

        # Lists of relative neighbour address (i.e. neighbours of (0,0)) which are used to calculate actual neighbours
        # through addition of coordinate values
        self.moore_relative = []
        self.von_neumann_relative = []
        self.populate_neighbour_tables(list_addresses)

    def populate_neighbour_tables(self, list_addresses):
        """
        Populates the full neighbour tables. First calculates the neighbour values relative to (0,0), then uses these
        against the actual addresses to create the lists
        e.g. von_neumann_relative[1] = [(-1,0), (0,-1), (0,1), (1,0)]
          then for, say, address (1,1) each value is added to (1,1) in turn to get
          [(0,1), (1,0), (1,2), (2,1)]
        :param list_addresses:
        :return:
        """
        # Initialise empty dictionaries
        self.moore_relative = dict()
        self.von_neumann_relative = dict()
        # Add an entry for each depth
        for d in range(1, self.max_depth + 1):
            self.von_neumann_relative[d] = []

        for depth in range(1, self.max_depth+1):
            # Get truth table values (e.g. depth 2 gives [-2,-1,0,1,2] for range_)
            range_ = range(-depth, depth + 1)
            # Use product to find all combinations for given depth and number of dimensions
            row = list(itertools.product(range_, repeat=self.dimensions))
            # Remove the 0 entry e.g. (0,0) for 2 dimensions
            row.remove((0,) * self.dimensions)

            reduced_row_moore = []
            self.von_neumann_relative[depth] = []
            for neighbour in row:
                # Calculate Manhattan distance and add to appropriate von Neumann table row
                manhattan_distance = int(sum([math.fabs(x) for x in neighbour]))
                if manhattan_distance <= self.max_depth and neighbour not in \
                        self.von_neumann_relative[manhattan_distance]:
                    self.von_neumann_relative[manhattan_distance].append(neighbour)
                # Check if one of coordinates = depth, if so then use for moore at this depth
                for x in neighbour:
                    if int(math.fabs(x)) == depth:
                        reduced_row_moore.append(neighbour)
                        break
            self.moore_relative[depth] = reduced_row_moore

        # Having calculated the relative Moore and VN addresses, apply these to every address to get a list of neighbour
        # addresses at varying depths for each address
        for address in list_addresses:
            self.moore_neighbours[address] = dict()
            self.von_neumann_neighbours[address] = dict()
            # Loop through each depth
            for depth in range(1, self.max_depth+1):
                # Apply the Moore neighbours and write to dictionary
                table = self.moore_relative[depth]
                output = [tuple(address[j] + table[i][j] for j in range(len(address))) for i in range(len(table))]
                self.moore_neighbours[address][depth] = output
                # Apply the Von neumann neighbours and write to dictionary
                table = self.von_neumann_relative[depth]
                output = [tuple(address[j] + table[i][j] for j in range(len(address))) for i in range(len(table))]
                self.von_neumann_neighbours[address][depth] = output

    def configure_neighbourhood_for_halo(self, halo_addresses):
        """
        Apply the neighbour calculation to the addresses in the halo (needed for calculation of caseum in neighbourhoods
        to reduce diffusion
        :param halo_addresses:
        :return:
        """
        # Loop through each halo address
        for address in halo_addresses:
            # Initially empty lists
            self.moore_neighbours[address] = dict()
            self.von_neumann_neighbours[address] = dict()
            # Loop through each depth up to max depth
            for depth in range(1, self.max_depth + 1):
                # Apply the Moore neighbours and write to dictionary
                table = self.moore_relative[depth]
                output = [tuple(address[j] + table[i][j] for j in range(len(address))) for i in range(len(table))]
                self.moore_neighbours[address][depth] = output
                # Apply the VN neighbours and write to dictionary
                table = self.von_neumann_relative[depth]
                output = [tuple(address[j] + table[i][j] for j in range(len(address))) for i in range(len(table))]
                self.von_neumann_neighbours[address][depth] = output


# TODO - COMP - not too keen on having this as a separate class, could just as easily be in the automaton class
class EventHandler:
    """
    Receives events (which have been deemed acceptable to perform) and performs them, altering agent values, adding and
    deleting where appropriate. The agents will be then be persisted onto the work grid.
    Because events can cross boundaries, need to check where each event affects. So check that the addresses affected
    lie on the grid.
    e.g. A macrophage moving from (a) to (b). (a) is on tile 0, (b) is on tile 1. When processed, the event is split
    into 2 events - the removal of macrophage from (a) on tile 0 and the addition of macrophage to (b) on tile 1. So
    event handler needs to ensure it only completes the actions of the event that relate to its tile (by checking
    self.address_is_on_grid where necessary)
    Also, where possible, pull agents from the grid not the event - passing events back and forth may change where
    they're stored in memory. So if an agent moves from (a) to (b) use the reference that exists in cell (a) to locate
    the agent
    Updates to agents are made to the agents themselves and the lists of agents, not to the grid. When a time-step
    finishes, all agents are written onto the grid for the next time-step
    """

    def __init__(self):
        pass

    def handle_event(self, event):
        """
        Given an event, call the appropriate function to make changes to the grid
        :param event:
        :return:
        """
        # Call relevant function based on the type of event
        if isinstance(event, BacteriumReplication):
            self.process_bacterium_replication(event)
        elif isinstance(event, RecruitTCell):
            self.process_t_cell_recruitment(event)
        elif isinstance(event, RecruitMacrophage):
            self.process_macrophage_recruitment(event)
        elif isinstance(event, ChemoKillBacterium):
            self.process_chemo_kill_bacterium(event)
        elif isinstance(event, ChemoKillMacrophage):
            self.process_chemo_kill_macrophage(event)
        elif isinstance(event, TCellDeath):
            self.process_t_cell_death(event)
        elif isinstance(event, TCellMovement):
            self.process_t_cell_movement(event)
        elif isinstance(event, TCellKillsMacrophage):
            self.process_t_cell_kill_macrophage(event)
        elif isinstance(event, MacrophageDeath):
            self.process_macrophage_death(event)
        elif isinstance(event, MacrophageMovement):
            self.process_macrophage_movement(event)
        elif isinstance(event, MacrophageKillsBacterium):
            self.process_macrophage_kills_bacterium(event)
        elif isinstance(event, MacrophageChangesState):
            self.process_macrophage_state_change(event)
        elif isinstance(event, BacteriumStateChange):
            self.process_bacterium_state_change(event)
        elif isinstance(event, MacrophageBursting):
            self.process_macrophage_bursting(event)
        else:
            raise Exception('Event ', type(event), 'not handled')

    def process_bacterium_replication(self, event):
        """
        Add a new bacteria to the grid
        :param event:
        :return:
        """
        # Is the new bacterium address on the grid
        if self.address_is_on_grid(event.new_bacterium_address):
            # Add new bacterium to grid
            self.add_bacterium(event.new_bacterium_address, event.new_metabolism)
        # Is the original bacteria address on the grid
        if self.address_is_on_grid(event.original_bacterium_address):
            # Pull the bacterium reference from the grid
            bacterium = self.grid[event.original_bacterium_address]['contents']
            # Reset age
            bacterium.age = 0.0
            # Swap the division neighbourhood between Moore and von Neumann
            if bacterium.division_neighbourhood == 'mo':
                bacterium.division_neighbourhood = 'vn'
            else:
                bacterium.division_neighbourhood = 'mo'

    def process_t_cell_recruitment(self, event):
        """
        Recruit a new t-cell from blood vessel
        :param event:
        :return:
        """
        # Only process if address is on the grid
        if self.address_is_on_grid(event.t_cell_address):
            self.add_t_cell(event.t_cell_address)

    def process_macrophage_recruitment(self, event):
        """
        Recruit a new macrophage from blood vessel
        :param event:
        :return:
        """
        # Only process if address is on the grid
        if self.address_is_on_grid(event.macrophage_address):
            self.add_macrophage(event.macrophage_address, "resting")

    def process_chemo_kill_bacterium(self, event):
        """
        Chemotherapy kills a bacterium
        :param event:
        :return:
        """
        # Pull bacteria reference from grid and remove it from the list
        bacterium = self.grid[event.address]['contents']
        self.bacteria.remove(bacterium)

    def process_chemo_kill_macrophage(self, event):
        """
        Chemotherapy kills an infected macrophage
        :param event:
        :return:
        """
        # Pull the macrophage reference from the grid and remove it
        macrophage = self.grid[event.dependant_addresses[0]]['contents']
        self.macrophages.remove(macrophage)
        # Add caseum at the address to the caseum list
        self.caseum.append(event.dependant_addresses[0])

    def process_t_cell_death(self, event):
        """
        A t-cell dies (through age)
        :param event:
        :return:
        """
        # Pull t-cell reference from grid and remove from list
        t_cell_to_die = self.grid[event.address]['contents']
        self.t_cells.remove(t_cell_to_die)

    def process_t_cell_movement(self, event):
        """
        A t-cell moves to a new cell
        :param event:
        :return:
        """
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        # T-cell moving between 2 cells in the same tile
        if event.internal:
            t_cell = self.grid[from_address]['contents']
            t_cell.address = to_address
        elif self.address_is_on_grid(from_address):  # T-cell is moving to a new tile
            t_cell = self.grid[from_address]['contents']
            self.t_cells.remove(t_cell)
        elif self.address_is_on_grid(to_address):  # T-cell has arrived from another tile
            event.t_cell_to_move.address = to_address
            self.t_cells.append(event.t_cell_to_move)

    def process_t_cell_kill_macrophage(self, event):
        """
        T-cell moves to a cell and kills an infected macrophage (also destroying itself)
        :param event:
        :return:
        """
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        # Is the macrophage on the grid
        if self.address_is_on_grid(to_address):
            # Remove macrophage, add caseum to cell
            macrophage = self.grid[to_address]['contents']
            self.macrophages.remove(macrophage)
            self.caseum.append(to_address)
        # Is the t-cell on the grid
        if self.address_is_on_grid(from_address):
            # Remove t-cell
            t_cell = self.grid[from_address]['contents']
            self.t_cells.remove(t_cell)

    def process_macrophage_death(self, event):
        """
        Macrophage dies, through age
        :param event:
        :return:
        """
        # Resting or active die, infected/chronically infected turn to caseum
        macrophage_to_die = self.grid[event.address]['contents']
        if macrophage_to_die.state == 'infected' or macrophage_to_die.state == 'chronically_infected':
            self.caseum.append(macrophage_to_die.address)
        # Remove macrophage
        self.macrophages.remove(macrophage_to_die)

    def process_macrophage_movement(self, event):
        """
        Macrophage moves from one cell to an empty cell
        :param event:
        :return:
        """
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        # Macrophage moving between 2 cells in the same tile
        if event.internal:
            macrophage = self.grid[from_address]['contents']
            macrophage.address = to_address
        elif self.address_is_on_grid(from_address):  # Macrophage is moving to a new tile
            # Remove macrophage
            macrophage = self.grid[from_address]['contents']
            self.macrophages.remove(macrophage)
        elif self.address_is_on_grid(to_address):  # Macrophage has arrived from another tile
            # Add macrophage
            event.macrophage_to_move.address = to_address
            self.macrophages.append(event.macrophage_to_move)

    def process_macrophage_kills_bacterium(self, event):
        """
        Macrophage moves to a cell with a bacterium, ingests said bacterium
        :param event:
        :return:
        """

        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]

        # Macrophage moving between 2 cells in the same tile
        if event.internal:
            # Move macrophage
            macrophage = self.grid[from_address]['contents']
            macrophage.address = to_address
            # Remove bacterium
            bacterium = self.grid[to_address]['contents']
            self.bacteria.remove(bacterium)
            # Active macrophages phagocytose bacteria, other states increment their intracellular bacteria count
            if macrophage.state == 'resting' or macrophage.state == 'infected' or macrophage.state == \
                    'chronically_infected':
                event.macrophage_to_move.intracellular_bacteria += 1
        elif self.address_is_on_grid(from_address):  # Macrophage is moving to a new tile
            macrophage = self.grid[from_address]['contents']
            self.macrophages.remove(macrophage)
        elif self.address_is_on_grid(to_address):  # Macrophage has arrived from another tile
            event.macrophage_to_move.address = to_address
            self.macrophages.append(event.macrophage_to_move)
            bacterium = self.grid[to_address]['contents']
            self.bacteria.remove(bacterium)
            # Active macrophages phagocytose bacteria, other states increment their intracellular bacteria count
            if event.macrophage_to_move.state == 'resting' or event.macrophage_to_move.state == 'infected' or \
                    event.macrophage_to_move.state == 'chronically_infected':
                event.macrophage_to_move.intracellular_bacteria += 1

    def process_macrophage_state_change(self, event):
        """
        Macrophage changes its state (e.g. resting -> active)
        :param event:
        :return:
        """
        # Pulling the macrophage from the grid
        macrophage = self.grid[event.address]['contents']
        # Change the state
        macrophage.state = event.new_state

    def process_bacterium_state_change(self, event):
        """
        A bacterium changes its state (fast <-> slow, resting -> true/false)
        :param event:
        :return:
        """
        # Pull bacterium from grid
        bacterium = self.grid[event.address]['contents']
        # Change the relevant state
        if event.type_of_change == 'metabolism':
            bacterium.metabolism = event.new_value
        elif event.type_of_change == 'resting':
            bacterium.resting = event.new_value

    def process_macrophage_bursting(self, event):
        """
        A macrophage bursts, spreading bacteria into the local neighbourhood
        :param event:
        :return:
        """
        # Is the macrophage on the grid
        if self.address_is_on_grid(event.macrophage_address):
            # Pull the macrophage from the grid
            macrophage_to_burst = self.grid[event.macrophage_address]['contents']
            # Remove macrophage from list and add caseum in its place
            self.macrophages.remove(macrophage_to_burst)
            self.caseum.append(macrophage_to_burst.address)

        # For each bacteria address, check it's an allowed address and the address is on this tile
        for i in event.bacteria_addresses:
            if i in event.impacted_addresses_allowed and self.address_is_on_grid(i):
                # Add a new slow bacteria there
                self.add_bacterium(i, 'slow')


class Automaton(Tile, Neighbourhood, EventHandler):

    def __init__(self, shape, tile_id, attributes, parameters, blood_vessels, fast_bacteria=None, slow_bacteria=None,
                 macrophages=None):
        Tile.__init__(self, shape, attributes)
        Neighbourhood.__init__(self, len(shape), parameters['max_depth'], self.list_grid_addresses)
        EventHandler.__init__(self)

        # Max depth must be +1 or more greater than the caseum distance as we need to work out diffusion rates one
        # cell deep into the halo
        assert parameters['caseum_distance_to_reduce_diffusion'] + 1 <= parameters['max_depth']

        self.tile_id = tile_id
        self.parameters = parameters
        self.time = 0

        # Local maxima (set internally)
        self.max_oxygen_local = 0.0
        self.max_chemotherapy_local = 0.0
        self.max_chemokine_local = 0.0
        # Global maxima (set externally)
        self.max_oxygen_global = 0.0
        self.max_chemotherapy_global = 0.0
        self.max_chemokine_global = 0.0
        self.number_of_bacteria_global = 0

        self.blood_vessels = []
        self.agents = []
        self.bacteria = []
        self.macrophages = []
        self.t_cells = []
        self.potential_events = []
        self.caseum = []

        # Initialise blood vessels and O2 levels
        self.initialise_blood_vessels(blood_vessels)
        self.initialise_oxygen_levels()

        # Calculate chemo schedule
        self.chemo_schedule1_start = np.random.randint(self.parameters['chemotherapy_schedule1_start_lower'],
                                                       self.parameters['chemotherapy_schedule1_start_upper'])

        # Copy grid to work grid
        self.create_work_grid()

        # Initialise agents onto the work grid
        self.initialise_bacteria(fast_bacteria, slow_bacteria)
        self.initialise_macrophages(macrophages)

        # Set up output file paths
        self.totalcell_test_file_path = str(self.tile_id) + '_totalcell_test.txt'
        self.contents_file_path = str(self.tile_id) + '_data_test.txt'
        self.oxygen_file_path = str(self.tile_id) + '_oxygen_test.txt'
        self.chemotherapy_file_path = str(self.tile_id) + '_chemo1.txt'
        self.chemokine_file_path = str(self.tile_id) + '_ckine.txt'
        self.type1_file_path = str(self.tile_id) + '_Type1.txt'
        self.type1_r_file_path = str(self.tile_id) + '_Type1_R.txt'
        self.type2_file_path = str(self.tile_id) + '_Type2.txt'
        self.type2_r_file_path = str(self.tile_id) + '_Type2_R.txt'
        self.type3_file_path = str(self.tile_id) + '_Type3.txt'
        self.activemac_file_path = str(self.tile_id) + '_activemac.txt'
        self.restingmac_file_path = str(self.tile_id) + '_restingmac.txt'
        self.infectedmac_file_path = str(self.tile_id) + '_infectedmac.txt'
        self.chroninfectedmac_file_path = str(self.tile_id) + '_chroninfectedmac.txt'
        self.caseation_file_path = str(self.tile_id) + '_caseation.txt'
        self.total_file_path = str(self.tile_id) + '_Total.txt'
        self.intra_bac_file_path = str(self.tile_id) + '_intra_bac.txt'

        # Clear up any old output files
        # TODO - should be check file exists, if so remove
        try:
            os.remove(self.contents_file_path)
            os.remove(self.oxygen_file_path)
            os.remove(self.chemotherapy_file_path)
            os.remove(self.chemokine_file_path)
            os.remove(self.type1_file_path)
            os.remove(self.type1_r_file_path)
            os.remove(self.type2_file_path)
            os.remove(self.type2_r_file_path)
            os.remove(self.type3_file_path)
            os.remove(self.activemac_file_path)
            os.remove(self.restingmac_file_path)
            os.remove(self.infectedmac_file_path)
            os.remove(self.chroninfectedmac_file_path)
            os.remove(self.caseation_file_path)
            os.remove(self.total_file_path)
            os.remove(self.intra_bac_file_path)
        except OSError:
            pass

        self.type1_file = open(self.type1_file_path, 'a')
        self.type1_r_file = open(self.type1_r_file_path, 'a')
        self.type2_file = open(self.type2_file_path, 'a')
        self.type2_r_file = open(self.type2_r_file_path, 'a')
        self.type3_file = open(self.type3_file_path, 'a')
        self.activemac_file = open(self.activemac_file_path, 'a')
        self.restingmac_file = open(self.restingmac_file_path, 'a')
        self.infectedmac_file = open(self.infectedmac_file_path, 'a')
        self.chroninfectedmac_file = open(self.chroninfectedmac_file_path, 'a')
        self.caseation_file = open(self.caseation_file_path, 'a')
        self.total_file = open(self.total_file_path, 'a')
        self.total_cell_test_file = open(self.totalcell_test_file_path, 'a')
        self.intra_bac_file = open(self.intra_bac_file_path, 'a')

        self.contents_file = open(self.contents_file_path, 'a')
        self.oxygen_file = open(self.oxygen_file_path, 'a')
        self.chemotherapy_file = open(self.chemotherapy_file_path, 'a')
        self.chemokine_file = open(self.chemokine_file_path, 'a')

        # Swap the working grid with actual grid to start process
        self.swap_grids()

    def close_files(self):
        self.type1_file.close()
        self.type1_r_file.close()
        self.type2_file.close()
        self.type2_r_file.close()
        self.type3_file.close()
        self.activemac_file.close()
        self.restingmac_file.close()
        self.infectedmac_file.close()
        self.chroninfectedmac_file.close()
        self.caseation_file.close()
        self.total_file.close()
        self.total_cell_test_file.close()
        self.intra_bac_file.close()

        self.contents_file.close()
        self.oxygen_file.close()
        self.chemotherapy_file.close()
        self.chemokine_file.close()

    def initialise_blood_vessels(self, addresses):
        """
        Set the initial blood vessels. These will not change.
        :param addresses:
        :return:
        """
        for address in addresses:
            # Set the value on grid and add to blood vessel list
            self.set_attribute_grid(address, 'blood_vessel', self.parameters['blood_vessel_value'])
            self.blood_vessels.append(address)

    def initialise_oxygen_levels(self):
        """
        Set the initial oxygen levels at each blood vessel location
        :return:
        """
        for address in self.blood_vessels:
            self.set_attribute_grid(address, 'oxygen',
                                    self.parameters['initial_oxygen']*self.grid[address]['blood_vessel'])
            self.max_oxygen_local = max(self.max_oxygen_local, self.parameters['initial_oxygen'])

    def initialise_bacteria(self, fast_bacteria, slow_bacteria):
        """
        Set the initial bacteria
        :param fast_bacteria:
        :param slow_bacteria:
        :return:
        """
        for address in fast_bacteria:
            self.add_bacterium(address, 'fast')
        for address in slow_bacteria:
            self.add_bacterium(address, 'slow')

    def initialise_macrophages(self, addresses):
        """
        Set the initial (resting) macrophages
        :param addresses:
        :return:
        """
        for address in addresses:
            self.add_macrophage(address, 'resting')

    def update(self):
        """
        Performs the continuous update (diffusion of oxygen, chemokine, chemotherapy) and applies to the work grid,
        and creates a list of potential events to be executed.
        Events are based on the state of the current GRID (i.e. are not affected by the diffusion changes)
        :return:
        """
        # Increment time
        self.time += 1

        # ----------------------------
        # CONTINUOUS (Diffusion)
        # ----------------------------

        # Pre-processing (calculating diffusion rates)
        self.diffusion_pre_process()
        # In chemo window?
        # TODO - COMP - chemo should be a duration
        chemo = (self.chemo_schedule1_start / self.parameters['time_step']) <= self.time < \
            (self.parameters['chemotherapy_schedule1_end'] / self.parameters['time_step']) or \
            self.parameters['chemotherapy_schedule2_start'] / self.parameters['time_step'] <= self.time

        # Reset local maximum values
        self.max_oxygen_local = 0.0
        self.max_chemotherapy_local = 0.0
        self.max_chemokine_local = 0.0

        # Loop through all cells
        self.diffusion(chemo)

        # ----------------------------
        # DISCRETE (Agents)
        # ----------------------------

        # Performs the necessary pre-condition checks for each event type and updates the list of potential events to
        # take place this time-step

        # Reset list
        self.potential_events = []

        self.bacteria_replication()

        self.t_cell_recruitment()

        self.macrophage_recruitment()

        self.chemotherapy_killing_bacteria()

        self.chemotherapy_killing_macrophages()

        self.t_cell_processes()

        self.macrophage_processes()

        self.macrophage_state_changes()

        self.bacteria_state_changes()

        self.reorder_events()

    def get_total_bacteria(self):
        return len(self.bacteria) + sum([m.intracellular_bacteria for m in self.macrophages])

    def diffusion_pre_process(self):
        """
        Diffusion of oxygen and chemo reduces inside a granuloma. Checks each cells neighbourhood for caseum, and if the
        total is above a threshold, oxygen and chemo diffusion are reduced accordingly.
        :return:
        """

        # Using the caseum list, record which cells are affected (i.e. within the neighbourhood). If a cell is affected
        # by enough caseum cells (over the threshold), reduce it's diffusion
        affected_addresses = []
        # Make a copy of the caseum addresses (don't directly reference as it will be added to with halo caseum)
        caseum_addresses = list(self.caseum)

        # Find caseum in the halo, add to the caseum list
        for address in self.list_halo_addresses:
            if self.grid[address] is not None and self.grid[address]['contents'] == 'caseum':
                caseum_addresses.append(address)

        # Process all caseum addresses
        for address in caseum_addresses:
            # Check the nighbourhood of the cell, up to the pre-set depth
            for depth in range(1, int(self.parameters['caseum_distance_to_reduce_diffusion']+1)):
                # Record each affected neighbour in the list (can be duplicates in list)
                neighbours = self.moore_neighbours[address][depth]
                for neighbour in neighbours:
                    affected_addresses.append(neighbour)

        # Count how many times each address appears in the affected list
        counted = Counter(affected_addresses)

        # TODO - COMP - is looping through every address inefficient? Once a cell has reduced diffusion it never reverts
        #  back. Also, we know which cells we're going to process from counted.

        # Loop through every address
        for address in self.list_grid_addresses:
            grid_cell = self.grid[address]
            # Get initial diffusion rates
            oxygen_diffusion = self.parameters['oxygen_diffusion']
            chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']
            # If the address has been affected and the number of affectations exceeds the threshold
            if address in counted and counted[address] >= \
                    self.parameters['caseum_threshold_to_reduce_diffusion']:
                # Reduce the diffusion rates at the cell
                oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                # Reduce the oxygen from source value
                if grid_cell['blood_vessel'] > 0.0:
                    grid_cell['blood_vessel'] = self.parameters['blood_vessel_value'] / \
                                                self.parameters['oxygen_diffusion_caseum_reduction']

            # Need to set the values on the current grid
            grid_cell['oxygen_diffusion_rate'] = oxygen_diffusion
            grid_cell['chemotherapy_diffusion_rate'] = chemotherapy_diffusion

        # Repeat the above, but for the halo addresses
        for halo_address in self.halo_depth1:
            halo_cell = self.grid[halo_address]
            if halo_cell is not None:

                # Get initial rates
                oxygen_diffusion = self.parameters['oxygen_diffusion']
                chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

                if halo_address in counted and counted[halo_address] >= \
                        self.parameters['caseum_threshold_to_reduce_diffusion']:
                    oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                    chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                    # Reduce the oxygen from source value
                    if halo_cell['blood_vessel'] > 0.0:
                        halo_cell['blood_vessel'] /= self.parameters['oxygen_diffusion_caseum_reduction']

                # Need to set the values on the halo
                halo_cell['oxygen_diffusion_rate'] = oxygen_diffusion
                halo_cell['chemotherapy_diffusion_rate'] = chemotherapy_diffusion

    def diffusion(self, chemo):
        """
        Oxygen, chemotherapy and chemokine diffuse through cells. Calculate new rates from the current values and set
        new values on the work grid.
        :param chemo:
        :return:
        """

        # Loop through every address, compare the value in the immediate von Neumann neighbourhood to calculate
        # diffusion. Then amend values based on agents/blood vessels
        for address in self.list_grid_addresses:

            work_grid_cell = self.work_grid[address]

            # Clear the work grid - agents will be added as they are processed
            work_grid_cell['contents'] = 0.0

            cell = self.grid[address]
            neighbour_addresses = self.von_neumann_neighbours[address][1]
            # neighbours = [self.grid[neighbour_address] for neighbour_address in neighbour_addresses]
            above, left, right, below = [self.grid[neighbour_address] for neighbour_address in neighbour_addresses]

            # Initialise expression
            oxygen_expression = 0
            chemotherapy_expression = 0
            chemokine_expression = 0

            # # Diffusion into neighbours
            # for neighbour in neighbours:
            #     # Only process if not boundary
            #     if neighbour is not None:
            #         oxygen_neighbour_diffusion = neighbour['oxygen_diffusion_rate']
            #         oxygen_expression += ((cell['oxygen_diffusion_rate'] + oxygen_neighbour_diffusion) / 2 * (
            #             neighbour['oxygen'] - cell['oxygen'])) / (
            #                              self.parameters['spatial_step']**2)
            #         if chemo:
            #             chemotherapy_neighbour_diffusion = neighbour['chemotherapy_diffusion_rate']
            #             chemotherapy_expression += ((cell['chemotherapy_diffusion_rate'] +
            #                 chemotherapy_neighbour_diffusion)
            #                 / 2 * (neighbour['chemotherapy'] - cell['chemotherapy'])) / \
            #                 (self.parameters['spatial_step']**2)
            #
            #         chemokine_neighbour_diffusion = self.parameters['chemokine_diffusion']
            #         chemokine_expression += ((self.parameters['chemokine_diffusion'] + chemokine_neighbour_diffusion)
            #             / 2 * (neighbour['chemokine'] - cell['chemokine'])) / (
            #             self.parameters['spatial_step'] * self.parameters['spatial_step'])

            if not None in [above, left, right, below]:
                oxygen_expression += ((((cell['oxygen_diffusion_rate'] + below['oxygen_diffusion_rate'])/2) *
                                        (below['oxygen'] - cell['oxygen'])
                                        - ((cell['oxygen_diffusion_rate'] + above['oxygen_diffusion_rate'])/2) *
                                        (cell['oxygen'] - above['oxygen']))
                                        / self.parameters['spatial_step']**2) + \
                                        ((((cell['oxygen_diffusion_rate'] + right['oxygen_diffusion_rate'])/2) *
                                        (right['oxygen'] - cell['oxygen'])
                                        - ((cell['oxygen_diffusion_rate'] + left['oxygen_diffusion_rate'])/2) *
                                        (cell['oxygen'] - left['oxygen']))
                                        / self.parameters['spatial_step']**2)
            # TOP LEFT
            elif above is None and left is None and not None in [right, below]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                right['oxygen'] - 2 * cell['oxygen'] + right['oxygen']) / (self.parameters['spatial_step']**2)
            # TOP
            elif above is None and not None in [left, right, below]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                left['oxygen'] - 2 * cell['oxygen'] + right['oxygen']) / (self.parameters['spatial_step']**2)
            # TOP RIGHT
            elif above is None and right is None and not None in [left, below]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + below['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                left['oxygen'] - 2 * cell['oxygen'] + left['oxygen']) / (self.parameters['spatial_step']**2)
            # LEFT
            elif left is None and not None in [above, right, below]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + below['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                right['oxygen'] - 2 * cell['oxygen'] + right['oxygen']) / (self.parameters['spatial_step']**2)
            # RIGHT
            elif right is None and not None in [above, left, below]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (below['oxygen'] - 2 * cell['oxygen'] + above['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                left['oxygen'] - 2 * cell['oxygen'] + left['oxygen']) / (self.parameters['spatial_step']**2)
            # BOTTOM LEFT
            elif left is None and below is None and not None in [above, right]:
                oxygen_expression += cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                right['oxygen'] - 2 * cell['oxygen'] + right['oxygen']) / (self.parameters['spatial_step']**2)
            # BOTTOM
            elif below is None and not None in [above, left, right]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                left['oxygen'] - 2 * cell['oxygen'] + right['oxygen']) / (self.parameters['spatial_step']**2)
            # BOTTOM RIGHT
            elif right is None and below is None and not None in [above, left]:
                oxygen_expression  += cell['oxygen_diffusion_rate'] * (above['oxygen'] - 2 * cell['oxygen'] + above['oxygen']) / (self.parameters['spatial_step']**2) + cell['oxygen_diffusion_rate'] * (
                left['oxygen'] - 2 * cell['oxygen'] + left['oxygen']) / (self.parameters['spatial_step']**2)



            # Amendments based on cell contents

            # OXYGEN

            # Add oxygen entering through blood vessel
            oxygen_expression += (self.parameters['oxygen_from_source'] * cell['blood_vessel'])
            # If there is bacteria in cell, then oxygen is taken up by bacteria so remove
            if cell['contents'] != 0.0 and isinstance(cell['contents'], Bacterium):
                oxygen_expression -= self.parameters['oxygen_uptake_from_bacteria'] * cell['oxygen']
            # Calculate new level
            new_oxygen = cell['oxygen'] + self.parameters['time_step'] * oxygen_expression
            work_grid_cell['oxygen'] = new_oxygen
            # Overwrite the maximum oxygen value if larger
            if new_oxygen > self.max_oxygen_local:
                self.max_oxygen_local = new_oxygen
            #self.max_oxygen_local = max(self.max_oxygen_local, new_oxygen)

            # CHEMOTHERAPY

            if chemo:
                # Release of chemotherapy from blood vessel
                chemotherapy_expression += (self.parameters['chemotherapy_from_source'] * cell['blood_vessel'])
                # Chemotherapy decay
                chemotherapy_expression -= self.parameters['chemotherapy_decay'] * cell['chemotherapy']
                # Calculate new level
                new_chemotherapy = cell['chemotherapy'] + self.parameters['time_step'] * chemotherapy_expression
                if new_chemotherapy > self.max_chemotherapy_local:
                    self.max_chemotherapy_local = new_chemotherapy
                # self.max_chemotherapy_local = max(self.max_chemotherapy_local, new_chemotherapy)
            else:
                self.max_chemotherapy_local = 0.0
                new_chemotherapy = 0.0
            work_grid_cell['chemotherapy'] = new_chemotherapy

            # CHEMOKINE

            # Release of chemokine by bacteria
            if cell['contents'] != 0.0 and isinstance(cell['contents'], Bacterium):
                chemokine_expression += self.parameters['chemokine_from_bacteria']
            # Release of chemokine by (non-resting) macrophages
            if cell['contents'] != 0.0 and isinstance(cell['contents'], Macrophage) and \
                            cell['contents'].state != 'resting':
                chemokine_expression += self.parameters['chemokine_from_macrophage']
            # Chemokine decay
            chemokine_expression -= self.parameters['chemokine_decay'] * cell['chemokine']
            # Calculate new level
            new_chemokine = cell['chemokine'] + self.parameters['time_step'] * chemokine_expression
            if new_chemokine > self.max_chemokine_local:
                self.max_chemokine_local = new_chemokine
            # self.max_chemokine_local = max(self.max_chemokine_local, new_chemokine)
            work_grid_cell['chemokine'] = new_chemokine

    def bacteria_replication(self):
        """
        Bacteria replicate (produce a new bacterium agent) once they reach a certain age.
        :return:
        """
        # Loop through every bacteria, check age against a (stochastic) threshold, generate event if age is higher than
        # threshold
        for bacterium in self.bacteria:
            # Increment age
            bacterium.age += self.parameters['time_step']

            # Skip if the bacterium is resting
            if bacterium.resting:
                continue

            if bacterium.metabolism == 'fast':
                maximum = self.parameters['bacteria_replication_fast_upper']
                minimum = self.parameters['bacteria_replication_fast_lower']
            else:  # Slow
                maximum = self.parameters['bacteria_replication_slow_upper']
                minimum = self.parameters['bacteria_replication_slow_lower']

            replication_time = np.random.randint(minimum, maximum) / self.parameters['time_step']

            # TODO - MED - why is this a modulo and not a > ?
            # If the time is sufficient enough, bacteria can replicate
            if self.time % replication_time == 0:

                # Look for free neighbours
                free_neighbours = []
                # TODO - COMP - maybe 4 shouldn't be hard coded?
                for depth in range(1, 4):
                    # Pull the neighbours from the appropriate neighbourhood
                    if bacterium.division_neighbourhood == 'mo':
                        neighbours = self.moore_neighbours[bacterium.address][depth]
                    else:
                        neighbours = self.von_neumann_neighbours[bacterium.address][depth]
                    # Find a free neighbour (not a blood vessel and contents == 0.0)
                    for neighbour_address in neighbours:
                        neighbour = self.grid[neighbour_address]
                        if neighbour is not None and neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                            free_neighbours.append(neighbour_address)
                    # If a free neighbour found, don't look at greater depths
                    if len(free_neighbours) > 0:
                        break
                # A free neighbour has not been found anywhere
                if len(free_neighbours) == 0:
                    # Bacterium will change to resting state (quorum sensing)
                    new_event = BacteriumStateChange(bacterium.address, 'resting', True)
                    self.potential_events.append(new_event)
                else:  # Free space found
                    # Pick a free neighbour at random
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                    # Determine if the event crosses a tile boundary
                    internal = self.address_is_on_grid(neighbour_address)
                    # Create event and add to list of potential events
                    new_event = BacteriumReplication(bacterium.address, neighbour_address, bacterium.metabolism,
                                                     internal)
                    self.potential_events.append(new_event)

    def t_cell_recruitment(self):
        """
        Once bacteria over entire system reach a threshold, t-cells enter the system. Creates an event to add a t-cell
        to a cell next to a blood vessel
        :return:
        """
        # When global amount of bacteria exceeds threshold
        if self.number_of_bacteria_global >= self.parameters['bacteria_threshold_for_t_cells']:
            # Each blood vessel
            for bv_address in self.blood_vessels:
                # Generate event if probability according to parameters
                r = np.random.randint(1, 101)
                if r <= self.parameters['t_cell_recruitment_probability']:
                    # Get von Neumann neighbours
                    neighbours = self.von_neumann_neighbours[bv_address][1]
                    # Reduce neighbours to those which are empty and have high enough cheokine level
                    free_neighbours = []
                    for neighbour_address in neighbours:
                        neighbour = self.grid[neighbour_address]
                        if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0 \
                                and self.chemokine_scale(neighbour_address) > \
                                self.parameters['chemokine_scale_for_t_cell_recruitment']:
                            free_neighbours.append(neighbour_address)
                    # Check there is free space
                    if len(free_neighbours) > 0:
                        # Pick one of the neighbours
                        neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                        # Create event
                        internal = self.address_is_on_grid(neighbour_address)
                        new_event = RecruitTCell(neighbour_address, internal)
                        self.potential_events.append(new_event)

    def macrophage_recruitment(self):
        """
        Each step for each source vessel, there is a probability that macrophage will be recruited
        :return:
        """
        if self.number_of_bacteria_global >= self.parameters['bacteria_threshold_for_macrophage_recruitment']:
            chemokine_threshold = self.parameters['chemokine_scale_for_macrophage_recruitment_above_threshold']
        else:
            chemokine_threshold = self.parameters['chemokine_scale_for_macrophage_recruitment_below_threshold']

        # Loop through each blood vessel
        for bv_address in self.blood_vessels:
            # Generate event with probability based on parameters
            r = np.random.randint(1, 101)
            if r <= self.parameters['macrophage_recruitment_probability']:
                # Get neighbours, then reduce to those that are free and have sufficient chemokine scale
                neighbours = self.von_neumann_neighbours[bv_address][1]
                free_neighbours = []
                for neighbour_address in neighbours:
                    neighbour = self.grid[neighbour_address]
                    if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0 and \
                            self.chemokine_scale(neighbour_address) > chemokine_threshold:
                        free_neighbours.append(neighbour_address)

                if len(free_neighbours) > 0:
                    # Pick one of the neighbours
                    chosen_neighbour = free_neighbours[np.random.randint(len(free_neighbours))]
                    # Create event
                    internal = self.address_is_on_grid(chosen_neighbour)
                    new_event = RecruitMacrophage(chosen_neighbour, internal)
                    self.potential_events.append(new_event)

    def chemotherapy_killing_bacteria(self):
        """
        Chemotherapy destroys bacterium if the level is high enough
        :return:
        """
        # Loop through all bacteria
        for bacterium in self.bacteria:
            # Check chemotherapy scale against relevant parameter based on metabolism
            chemo_scale = self.chemotherapy_scale(bacterium.address)
            if (bacterium.metabolism == 'fast' and chemo_scale >
                self.parameters['chemotherapy_scale_for_kill_fast_bacteria']) \
                    or \
                    (bacterium.metabolism == 'slow' and chemo_scale >
                        self.parameters['chemotherapy_scale_for_kill_slow_bacteria']):
                # Scale is high enough, so create event to destroy bacterium
                new_event = ChemoKillBacterium(bacterium.address)
                self.potential_events.append(new_event)

    def chemotherapy_killing_macrophages(self):
        """
        Chemotherapy destroys infected macrophage if the level is high enough
        :return:
        """
        # Loop through all macrophages
        for m in self.macrophages:
            # Check scale against parameter threshold (only for infected or chronically infected macrophages)
            chemo_scale = self.chemotherapy_scale(m.address)
            if ((m.state == 'infected' or m.state == 'chronically_infected') and chemo_scale >
                    self.parameters['chemotherapy_scale_for_kill_macrophage']):
                # Scale if high enough so create event
                new_event = ChemoKillMacrophage(m.address)
                self.potential_events.append(new_event)

    def t_cell_processes(self):
        """
        T-cells movement, death and apoptosis of other agents
        :return:
        """
        # T-cells only move after set period of time
        if self.time % self.parameters['t_cell_movement_time'] == 0:

            # Loop through all T-cells
            for t_cell in self.t_cells:
                # Increment age
                t_cell.age += self.parameters['time_step']
                # Stochastic age threshold
                age_threshold = np.random.randint(0, self.parameters['t_cell_age_threshold'])
                # T-CELL DEATH
                # If age > threshold, t-cell dies
                if t_cell.age >= age_threshold:
                    new_event = TCellDeath(t_cell.address)
                    self.potential_events.append(new_event)
                else:  # T-CELL MOVE
                    # T-cells move in biased random walk. Determine if move will be random based on probability in
                    # parameters
                    random_move = False
                    prob_random_move = np.random.randint(1, 101)
                    if prob_random_move <= self.parameters['t_cell_random_move_probability']:
                        random_move = True
                    # Get neighbours
                    neighbours = self.moore_neighbours[t_cell.address][1]
                    # If a random move, pick a neighbour at random
                    if random_move:
                        # Remove neighbours not on system
                        possible_neighbours = [n for n in neighbours if self.grid[n] is not None]
                        index = np.random.randint(0, len(possible_neighbours))
                        chosen_neighbour_address = possible_neighbours[index]
                    else: # Pick the neighbour with the highest chemokine level
                        chosen_index = self.find_max_chemokine_neighbour(neighbours)[0]
                        chosen_neighbour_address = neighbours[chosen_index]

                    # Check if leaving the grid
                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    # Get neighbour
                    neighbour = self.grid[chosen_neighbour_address]
                    # Check neighbour is empty, then move T-cell there
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = TCellMovement(t_cell, t_cell.address, chosen_neighbour_address, internal)
                        self.potential_events.append(new_event)
                    # Else if the address contains an infected macrophage, then t-cell may kill it
                    elif isinstance(neighbour['contents'], Macrophage) and (neighbour['contents'].state == 'infected'
                                                                            or neighbour[
                            'contents'].state == 'chronically_infected'):
                        # T-cell killing based on parameter probability
                        prob_t_cell_killing = np.random.randint(1, 101)
                        if prob_t_cell_killing <= self.parameters['t_cell_kills_macrophage_probability']:
                            new_event = TCellKillsMacrophage(t_cell, t_cell.address, chosen_neighbour_address, internal)
                            self.potential_events.append(new_event)

    def macrophage_processes(self):
        """
        Macrophages move, die and ingest bacteria
        :return:
        """
        # Loop through macrophages
        for macrophage in self.macrophages:
            # Increment age
            macrophage.age += self.parameters['time_step']
            # Different events/movement rates/death rates depending on state
            if macrophage.state == 'resting':
                # Death by age is stochastic
                random_macrophage_age = np.random.randint(0, self.parameters['resting_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age == 0:
                    # Create an event
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue
                # Within a set time for movement
                if self.time % self.parameters['resting_macrophage_movement_time'] == 0:
                    # Chemokine moves on random biased walk. Random move with probability based on parameters, if
                    # highest chemokine scale at neighbours does not exceed threshold, then also random move
                    neighbours = [n for n in self.moore_neighbours[macrophage.address][1] if self.grid[n] is not None]
                    chosen_index, max_chemokine_scale = self.find_max_chemokine_neighbour(neighbours)
                    # Generate random number for probability of random move
                    prob_random_move = np.random.randint(1, 101)
                    random_move = False
                    if prob_random_move <= self.parameters['prob_resting_macrophage_random_move'] \
                            or max_chemokine_scale <= \
                            self.parameters['minimum_chemokine_for_resting_macrophage_movement']:
                        random_move = True
                    # Pick the neighbour to move to, either random or highest chemokine scale
                    if random_move:
                        chosen_neighbour_address = neighbours[np.random.randint(0, len(neighbours))]
                    else:
                        chosen_neighbour_address = neighbours[chosen_index]
                    # Check if leaving the grid
                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    neighbour = self.grid[chosen_neighbour_address]
                    # If neighbour is empty, create a move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    # If neighbour contains a bacterium, ingest it
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)
            # Active macrophage processes
            elif macrophage.state == 'active':
                # Active macrophages die after a set time (not stochastic)
                if macrophage.age > self.parameters['active_macrophage_age_limit']:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue
                # Set time for macrophage movement
                if self.time % self.parameters['active_macrophage_movement_time'] == 0:
                    # Active macrophages always move to highest chemokine neighbour
                    neighbours = [n for n in self.moore_neighbours[macrophage.address][1] if self.grid[n] is not None]
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]
                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    neighbour = self.grid[chosen_neighbour_address]
                    # If cell to move to has a bacterium
                    if isinstance(neighbour['contents'], Bacterium):
                        # Macrophages ingests with set probability (active macrophages will destroy)
                        prob_macrophage_kill = np.random.randint(1, 101)
                        # Probabilites differ based on bacterium metabolism
                        if (neighbour['contents'].metabolism == 'fast' and prob_macrophage_kill <= self.parameters[
                                'prob_active_macrophage_kill_fast_bacteria']) or (
                                neighbour['contents'].metabolism == 'slow' and prob_macrophage_kill <=
                                self.parameters['prob_active_macrophage_kill_slow_bacteria']):
                            new_event = MacrophageKillsBacterium(macrophage, macrophage.address,
                                                                 chosen_neighbour_address, internal)
                            self.potential_events.append(new_event)
                    # Cell is empty so create a move event
                    elif neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
            # Infected Macrophage processes
            elif macrophage.state == 'infected':
                # Death is stochastic
                random_macrophage_age = np.random.randint(0, self.parameters['infected_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age == 0:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue
                # Move after certain time
                if self.time % self.parameters['infected_macrophage_movement_time'] == 0:
                    # Infected move to highest chemokine neighbour
                    neighbours = [n for n in self.moore_neighbours[macrophage.address][1] if self.grid[n] is not None]
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]
                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    neighbour = self.grid[chosen_neighbour_address]
                    # Neighbour is empty, so move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    # Neighbour has a bacterium, so kill event
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)
            # Chronically infected macrophage processes
            elif macrophage.state == 'chronically_infected':
                # Stochastic death
                random_macrophage_age = np.random.randint(0,
                                                          self.parameters['chronically_infected_macrophage_age_limit'])
                if macrophage.age >= random_macrophage_age == 0:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue
                # Movement at set times
                if self.time % self.parameters['chronically_infected_macrophage_movement_time'] == 0:
                    # Move to highest chemokine scale neighbour
                    neighbours = [n for n in self.moore_neighbours[macrophage.address][1] if self.grid[n] is not None]
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]
                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    neighbour = self.grid[chosen_neighbour_address]
                    # Neighbour is empty, so move event
                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    # Neighbour has bacterium, so kill event
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)

    def macrophage_state_changes(self):
        """
        Macrophages change state based on chemokine levels and intracellular bacteria
        :return:
        """
        # Loop through macrophages
        for macrophage in self.macrophages:
            if macrophage.state == 'resting':
                # Resting to active if scale exceeds threshold
                if self.chemokine_scale(macrophage.address) > \
                        self.parameters['chemokine_scale_for_macrophage_activation'] \
                        and macrophage.intracellular_bacteria == 0:
                    new_event = MacrophageChangesState(macrophage.address, 'active')
                    self.potential_events.append(new_event)
                # Resting to infected if there is one intracellular bacteria
                elif macrophage.intracellular_bacteria == 1:
                    new_event = MacrophageChangesState(macrophage.address, 'infected')
                    self.potential_events.append(new_event)
            elif macrophage.state == 'active':
                # Active to resting if scale is low enough
                if self.chemokine_scale(macrophage.address) < \
                        self.parameters['chemokine_scale_for_macrophage_deactivation']:
                    new_event = MacrophageChangesState(macrophage.address, 'resting')
                    self.potential_events.append(new_event)
            elif macrophage.state == 'infected':
                # Infected to Chronically Infected if intracellular bacteria exceeds threshold
                if macrophage.intracellular_bacteria > self.parameters['bacteria_to_turn_chronically_infected']:
                    new_event = MacrophageChangesState(macrophage.address, 'chronically_infected')
                    self.potential_events.append(new_event)
            elif macrophage.state == 'chronically_infected':
                # Macrophage bursts if intracellular bacteria exceed threshold
                if macrophage.intracellular_bacteria == self.parameters['bacteria_to_burst_macrophage']:
                    internal = self.address_is_on_grid(macrophage.address)
                    # Loop through all neighbours (up to depth 3) and try to find enough to distribute bacteria on to
                    bacteria_addresses = []
                    for depth in range(1, 4):
                        neighbours = self.moore_neighbours[macrophage.address][depth]
                        # Shuffle the neighbours so we don't give priority
                        np.random.shuffle(neighbours)
                        for n in neighbours:
                            # Find empty neighbours
                            neighbour = self.grid[n]
                            if neighbour is not None and neighbour['contents'] == 0.0 and \
                                            neighbour['blood_vessel'] == 0.0:
                                bacteria_addresses.append(n)
                                # If an address is off the tile then event will be external
                                if not self.address_is_on_grid(n):
                                    internal = False
                            # Limit reached - break here stops checking other neighbours at this depth, also need to
                            # stop searching further depths
                            if len(bacteria_addresses) == self.parameters['bacteria_to_burst_macrophage']:
                                break
                        # Limit reached earlier so don't check other depths
                        if len(bacteria_addresses) == self.parameters['bacteria_to_burst_macrophage']:
                            break
                    # Macrophage bursting event
                    new_event = MacrophageBursting(macrophage.address, bacteria_addresses, internal)
                    self.potential_events.append(new_event)

    def bacteria_state_changes(self):
        """
        Bacteria switch between metabolism based on oxygen and switch resting true to false based on space
        :return:
        """
        # Loop through bacteria
        for bacterium in self.bacteria:
            # Metabolism change only happens later in process (after 2 hours)
            if self.time > 2 / self.parameters['time_step']:
                # Check if state change - different scales based on metabolism
                if bacterium.metabolism == 'fast' and self.oxygen_scale(bacterium.address) <= self.parameters[
                        'oxygen_scale_for_metabolism_change_to_slow']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'slow')
                    self.potential_events.append(new_event)
                elif bacterium.metabolism == 'slow' and self.oxygen_scale(bacterium.address) > self.parameters[
                        'oxygen_scale_for_metabolism_change_to_fast']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'fast')
                    self.potential_events.append(new_event)
            # If bacteria is resting, check if there is now space in neighbourhood, if so, revert to non-resting
            if bacterium.resting:
                space_found = False
                for depth in range(1, 4):
                    # Get neighbours
                    neighbours = self.moore_neighbours[bacterium.address][depth]
                    for n in neighbours:
                        # Is neighbour empty?
                        neighbour = self.grid[n]
                        if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0:
                            new_event = BacteriumStateChange(bacterium.address, 'resting', False)
                            self.potential_events.append(new_event)
                            space_found = True
                            # Don't check other neighbours
                            break
                    # Space found so don't check further depths
                    if space_found:
                        break

    def reorder_events(self):
        """
        Order the set of potential events
        :return:
        """
        # TODO - COMP - other methods - currently just random
        np.random.shuffle(self.potential_events)

    def process_events(self, events):
        """
        Given a list of acceptable events, pass to the Event Handler to update the grid
        :param events:
        :return:
        """
        # Handle each event through the event handler
        for event in events:
            self.handle_event(event)

        # Ensure all agents are put onto the work grid
        self.persist_agents()
        # Swap work grid with main grid
        self.swap_grids()

        # Record state
        if self.time % self.parameters['interval_to_record_results'] == 0.0:
            self.record_state()
        # Record counts
        self.record_counts()

    def set_max_oxygen_global(self, max_oxygen):
        """
        Set global maximum oxygen level (from external source)
        :param max_oxygen:
        :return:
        """
        self.max_oxygen_global = max_oxygen

    def set_max_chemotherapy_global(self, max_chemotherapy):
        """
        Set global maximum chemotherapy level (from external source)
        :param max_chemotherapy:
        :return:
        """
        self.max_chemotherapy_global = max_chemotherapy

    def set_max_chemokine_global(self, max_chemokine):
        """
        Set global maximum chemokine level (from external source)
        :param max_chemokine:
        :return:
        """
        self.max_chemokine_global = max_chemokine

    def set_global_bacteria_number(self, number):
        """
        Set global bacteria numbers (from external source)
        :param number:
        :return:
        """
        self.number_of_bacteria_global = number

    def oxygen_scale(self, address):
        """
        Oxygen level at cell as % of maximum global oxygen level
        :param address:
        :return:
        """
        if self.max_oxygen_global == 0.0:
            return 0.0
        else:
            return (self.grid[address]['oxygen'] / self.max_oxygen_global) * 100

    def chemotherapy_scale(self, address):
        """
        Maximum chemotherapy level at cell as % of global maximum chemotherapy
        :param address:
        :return:
        """
        if self.max_chemotherapy_global == 0.0:
            return 0.0
        else:
            return (self.grid[address]['chemotherapy'] / self.max_chemotherapy_global) * 100

    def chemokine_scale(self, address):
        """
        Chemokine level at cell as % of global maximum chemokine level
        :param address:
        :return:
        """
        if self.max_chemokine_global == 0.0:
            return 0.0
        else:
            return (self.grid[address]['chemokine'] / self.max_chemokine_global) * 100.0

    def add_bacterium(self, address, metabolism):
        """
        Add a new bacterium of set metabolism
        :param address:
        :param metabolism:
        :return:
        """
        new_bacterium = Bacterium(address, metabolism)
        self.bacteria.append(new_bacterium)
        self.set_attribute_work_grid(address, 'contents', new_bacterium)

    def add_macrophage(self, address, state):
        """
        Add a macrophage of set state
        :param address:
        :param state:
        :return:
        """
        new_macrophage = Macrophage(address, state)
        self.macrophages.append(new_macrophage)
        self.set_attribute_work_grid(address, 'contents', new_macrophage)

    def add_t_cell(self, address):
        """
        Add a T-cell at address
        :param address:
        :return:
        """
        new_t_cell = TCell(address)
        self.t_cells.append(new_t_cell)
        self.set_attribute_work_grid(address, 'contents', new_t_cell)

    def persist_agents(self):
        """
        For every agent in the lists, assign them to the cell contents using their addess attribute
        Ensures any agents who don't move, etc. remain on the system
        :return:
        """
        for b in self.bacteria:
            self.set_attribute_work_grid(b.address, 'contents', b)
        for m in self.macrophages:
            self.set_attribute_work_grid(m.address, 'contents', m)
        for t in self.t_cells:
            self.set_attribute_work_grid(t.address, 'contents', t)
        for c_address in self.caseum:
            self.set_attribute_work_grid(c_address, 'contents', 'caseum')

    def find_max_chemokine_neighbour(self, neighbours):
        """
        Given neighbour addresses, find the neighbour which has the highest level of chemokine
        :param neighbours:
        :return:
        """
        max_chemokine_scale = 0
        highest_indices = []
        for index in range(len(neighbours)):
            if self.grid[neighbours[index]] is not None:
                chemokine_scale = self.chemokine_scale(neighbours[index])
                if chemokine_scale > max_chemokine_scale:
                    max_chemokine_scale = chemokine_scale
                    highest_indices = [index]
                elif chemokine_scale == max_chemokine_scale:
                    highest_indices.append(index)

        # Tie-breaking. If just one pick it, else pick any one index at random
        choice = np.random.randint(0, len(highest_indices))
        chosen_index = highest_indices[choice]

        return chosen_index, max_chemokine_scale

    def record_state(self):
        """
        Record the attributes of each cell. Each attribute has a separate file.
        :return:
        """

        for address in self.list_grid_addresses:
            cell = self.grid[address]

            contents_number = 0.0

            # Record contents
            # 0.0 - Empty
            # 1.0 - Fast bacteria //  1.25 - fast bacteria resting
            # 2.0 - slow bacteria //  2.25 - slow bacteria resting
            # 3.0 - T-cell
            # 4.0 - Resting macrophage  //  5.0 - active macrophage
            # 6.0 - infected macrophage //  7.0 - chronically infected macrophage
            # 100.0 - caseum

            if isinstance(cell['contents'], Bacterium):
                if cell['contents'].metabolism == 'fast':
                    contents_number = 1.0
                else:
                    contents_number = 2.0
                if cell['contents'].resting:
                    contents_number += 0.25
            elif isinstance(cell['contents'], Macrophage):
                contents_number = 4.0
                if cell['contents'].state == 'active':
                    contents_number += 1.0
                elif cell['contents'].state == 'infected':
                    contents_number += 2.0
                elif cell['contents'].state == 'chronically_infected':
                    contents_number += 3.0
            elif isinstance(cell['contents'], TCell):
                contents_number = 3.0
            elif cell['contents'] == 'caseum':
                contents_number = 100.0
            self.contents_file.write(str(contents_number) + '\n')

            # Record oxygen
            self.oxygen_file.write(str(self.oxygen_scale(address)) + '\n')

            # Record chemotherapy
            self.chemotherapy_file.write(str(self.chemotherapy_scale(address)) + '\n')

            # Record chemokine
            self.chemokine_file.write(str(self.chemokine_scale(address)) + '\n')

    def record_counts(self):
        """
        Record the counts of various agents
        :return:
        """

        # Files for counts
        bacteria_count = len(self.bacteria)
        # Write fast bacteria numbers to file
        type1_count = len([n for n in self.bacteria if n.metabolism == 'fast' and not n.resting])
        self.type1_file.write(str(type1_count) + '\n')

        # Write fast resting bacteria numbers to file
        type1_r_count = len([n for n in self.bacteria if n.metabolism == 'fast' and n.resting])
        self.type1_r_file.write(str(type1_r_count) + '\n')

        # Write slow bacteria numbers to file
        type2_count = len([n for n in self.bacteria if n.metabolism == 'slow' and not n.resting])
        self.type2_file.write(str(type2_count) + '\n')

        # Write slow resting bacteria numbers to file
        type2_r_count = len([n for n in self.bacteria if n.metabolism == 'slow' and n.resting])
        self.type2_r_file.write(str(type2_r_count) + '\n')

        # Write t-cell numbers to file
        t_cell_count = len(self.t_cells)
        self.type3_file.write(str(t_cell_count) + '\n')

        # Write macrophage numbers to file
        macrophage_count = len(self.macrophages)

        activemac_count = 0
        restingmac_count = 0
        infectedmac_count = 0
        chroninfectedmac_count = 0
        for m in self.macrophages:
            if m.state == 'active':
                activemac_count += 1
            elif m.state == 'resting':
                restingmac_count += 1
            elif m.state == 'infected':
                infectedmac_count += 1
            elif m.state == 'chroninfected':
                chroninfectedmac_count += 1
        self.activemac_file.write(str(activemac_count) + '\n')
        self.restingmac_file.write(str(restingmac_count) + '\n')
        self.infectedmac_file.write(str(infectedmac_count) + '\n')
        self.chroninfectedmac_file.write(str(chroninfectedmac_count) + '\n')

        # Write caseum numbers to file
        caseation_count = len(self.caseum)
        self.caseation_file.write(str(caseation_count) + '\n')

        # Write total agent numbers to file
        total_count = bacteria_count + macrophage_count + t_cell_count + caseation_count
        self.total_file.write(str(total_count) + '\n')

        if self.time % self.parameters['interval_to_record_results'] == 0.0:
            self.total_cell_test_file.write(str(total_count) + '\n')

        # Write intracellular bacteria numbers to file
        intra_bac_count = sum([m.intracellular_bacteria for m in self.macrophages])
        self.intra_bac_file.write(str(intra_bac_count) + '\n')

# ------------------------------
# AGENTS
# ------------------------------


class Agent:

    def __init__(self, address):
        self.address = address
        self.age = 0


class Bacterium(Agent):

    def __init__(self, address, metabolism):
        self.metabolism = metabolism
        self.division_neighbourhood = 'mo'
        self.resting = False
        Agent.__init__(self, address)


class Macrophage(Agent):

    def __init__(self, address, state):
        self.state = state
        self.intracellular_bacteria = 0
        Agent.__init__(self, address)


class TCell(Agent):

    def __init__(self, address):
        Agent.__init__(self, address)
        pass


# ------------------------------
# EVENTS
# ------------------------------

class Event(object):
    def __init__(self, dependant_addresses, impacted_addresses, internal):
        self.dependant_addresses = dependant_addresses
        self.impacted_addresses_potential = impacted_addresses
        self.impacted_addresses_allowed = []
        self.internal = internal

    def clone(self, new_addresses):
        raise NotImplementedError


class BacteriumReplication(Event):
    def __init__(self, original_bacterium_address, new_bacterium_address, metabolism, internal):
        Event.__init__(self, [original_bacterium_address, new_bacterium_address],
                       [original_bacterium_address, new_bacterium_address], internal)
        self.new_bacterium_address = new_bacterium_address
        self.new_metabolism = metabolism
        self.original_bacterium_address = original_bacterium_address

    def clone(self, new_addresses):

        return BacteriumReplication(new_addresses[0], new_addresses[1], self.new_metabolism,
                                    self.internal)


class RecruitTCell(Event):
    def __init__(self, address, internal):
        self.t_cell_address = address
        Event.__init__(self, [address], [address], internal)

    def clone(self, new_addresses):
        return RecruitTCell(new_addresses[0], self.internal)


class RecruitMacrophage(Event):
    def __init__(self, address, internal):
        self.macrophage_address = address
        Event.__init__(self, [address], [address], internal)

    def clone(self, new_addresses):
        return RecruitMacrophage(new_addresses[0], self.internal)


class ChemoKillBacterium(Event):
    def __init__(self, address):
        self.address = address
        # Chemo killing is always internal
        Event.__init__(self, [address], [address], True)


class ChemoKillMacrophage(Event):
    def __init__(self, address):
        self.macrophage_to_kill_address = address
        # Chemo killing is always internal
        Event.__init__(self, [address], [address], True)


class TCellDeath(Event):
    def __init__(self, address):
        self.address = address
        # T-cell death is always internal
        Event.__init__(self, [address], [address], True)


class TCellMovement(Event):
    def __init__(self, t_cell_to_move, from_address, to_address, internal):
        self.t_cell_to_move = t_cell_to_move
        self.to_address = to_address

        Event.__init__(self, [from_address, to_address], [from_address, to_address], internal)

    def clone(self, new_addresses):
        return TCellMovement(self.t_cell_to_move, new_addresses[0], new_addresses[1], self.internal)


class TCellKillsMacrophage(Event):
    def __init__(self, t_cell, t_cell_address, macrophage_address, internal):
        self.t_cell_to_move = t_cell
        self.macrophage_address = macrophage_address
        Event.__init__(self, [t_cell_address, macrophage_address], [t_cell_address, macrophage_address], internal)

    def clone(self, new_addresses):
        return TCellKillsMacrophage(self.t_cell_to_move, new_addresses[0], new_addresses[1], self.internal)


class MacrophageDeath(Event):
    def __init__(self, address):
        self.address = address
        # Macrophage death is always internal
        Event.__init__(self, [address], [address], True)


class MacrophageMovement(Event):
    def __init__(self, macrophage_to_move, from_address, to_address, internal):
        self.macrophage_to_move = macrophage_to_move
        self.new_address = to_address
        Event.__init__(self, [from_address, to_address], [from_address, to_address], internal)

    def clone(self, new_addresses):
        return MacrophageMovement(self.macrophage_to_move, new_addresses[0], new_addresses[1], self.internal)


class MacrophageKillsBacterium(Event):
    def __init__(self, macrophage_to_move, macrophage_address, bacterium_address, internal):
        self.macrophage_to_move = macrophage_to_move
        self.bacterium_address = bacterium_address
        Event.__init__(self, [macrophage_address, bacterium_address], [macrophage_address, bacterium_address], internal)

    def clone(self, new_addresses):
        return MacrophageKillsBacterium(self.macrophage_to_move, new_addresses[0], new_addresses[1], self.internal)


class MacrophageChangesState(Event):
    def __init__(self, address, new_state):
        self.address = address
        self.new_state = new_state
        Event.__init__(self, [address], [address], True)


class BacteriumStateChange(Event):
    def __init__(self, address, type_of_change, new_value):
        self.address = address
        self.type_of_change = type_of_change
        self.new_value = new_value
        Event.__init__(self, [address], [address], True)


class MacrophageBursting(Event):
    def __init__(self, macrophage_address, bacteria_addresses, internal):
        self.macrophage_address = macrophage_address
        self.bacteria_addresses = bacteria_addresses
        impacted_addresses = [macrophage_address]
        impacted_addresses += bacteria_addresses
        Event.__init__(self, [macrophage_address], impacted_addresses, internal)