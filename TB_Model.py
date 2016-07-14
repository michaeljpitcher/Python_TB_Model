import numpy as np
import math
import itertools
from collections import Counter

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
Location - an integer value for a cell (e.g. in a 2D grid of shape [5,5], location 1 == address [0,1],
           location 8 == address [1,3]
Halo - the cells required by a tile for updating which are not a part of the tile (belong to other tiles)
Danger Zone - the cells in a tile which will be required by other tiles (i.e. are part of other tiles' halos)
'''


class Topology:

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local,
                 fast_bacteria_local, slow_bacteria_local, macrophages_local):

        self.number_of_tiles = reduce(lambda x, y: x * y, tile_arrangement)
        self.total_shape = np.array(total_shape)
        self.tile_shape = self.total_shape / tile_arrangement
        self.tile_arrangement = tile_arrangement

        self.automata = []

        # Create automata
        for i in range(self.number_of_tiles):
            automaton = Automaton(self.tile_shape, i, attributes, parameters, blood_vessel_local[i],
                                  fast_bacteria_local[i], slow_bacteria_local[i], macrophages_local[i])
            self.automata.append(automaton)

        # Get a list of addresses (relative to tile) that are required by tile but outside of it
        # Will be the same for all tiles (as they're same shape) so just pull from one
        self.external_addresses_required = self.get_external_addresses_required(self.automata[0],
                                                                                parameters['max_depth'])

        # Halo of depth 1 - needed for calculating diffusion rates
        # TODO - COMP - is there a better way of doing this?
        depth1_addresses = self.get_external_addresses_required(self.automata[0], 1)

        # Set the halo address and halo of depth 1 address
        for automaton in self.automata:
            automaton.halo_addresses = self.external_addresses_required
            automaton.halo_depth1 = depth1_addresses

    def get_external_addresses_required(self, automaton, depth):
        """
        The addresses of cells within the overall neighbourhood which don't belong to the given automaton
        :return:
        """
        external_addresses = []
        # Loop through every cell
        for i in range(0, automaton.grid.size):
            # Pull the addresses of cells in neighbourhood of this cell
            address = automaton.location_to_address(i)
            neighbours = []
            for j in range(1, depth+1):
                neighbours += automaton.neighbours_moore(address, j)
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

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_global,
                 fast_bacteria_global, slow_bacteria_global, macrophages_global):
        # Check it's two dimensions
        assert len(total_shape) == 2
        self.number_of_tiles = reduce(lambda a, q: a * q, tile_arrangement)
        self.total_shape = np.array(total_shape)
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

        # Use required global addresses to create danger zones on tiles
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
        """
        Turn a global address into a tile ID and a local address
        :param global_address:
        :return:
        """
        if global_address is None:
            return [None, None]

        x, y = global_address
        tile_rows, tile_cols = self.tile_shape

        # Add the tile ID - x modulo num of rows in a tile * number of tiles in width of the grid
        output = [divmod(x, tile_rows)[0] * self.tile_arrangement[1] + divmod(y, tile_cols)[0]]
        # Add the Local coordinates
        output.append([divmod(x, tile_rows)[1], divmod(y, tile_cols)[1]])

        return output

    def local_to_global(self, tile_id, local_address):
        """
        Turn a tile ID and a local address into a global address
        :param tile_id:
        :param local_address:
        :return:
        """
        x, y = local_address
        origin = self.origins[tile_id]
        # Normalise the address before it's returned
        return self.normalise_address([origin[0] + x, origin[1] + y])

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

    def __init__(self, shape, attributes):
        self.shape = shape
        self.attributes = attributes
        self.size = reduce(lambda x, y: x * y, shape)
        self.grid = self.create_grid(attributes)

    def create_grid(self, attributes):
        """
        Create a grid where cells are dictionaries of given attributes
        :param attributes:
        :return:
        """

        cells = []
        for i in range(self.size):
            cell = dict()
            for att in attributes:
                cell[att] = 0.0
            cells.append(cell)
        # Turn flat list into the necessary dimensional array
        return np.array(cells).reshape(self.shape)

    def create_work_grid(self):
        """
        Create a working grid placeholder by cloning the current grid
        :return:
        """
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
        """
        Get the cell values of cells in the danger zone
        :return:
        """
        danger_zone = []
        for address in self.danger_zone_addresses:
            # Get cell from the work grid
            danger_zone.append(self.grid[tuple(address)])

        return danger_zone

    def set_addresses_for_halo(self, addresses):
        self.halo_addresses = addresses

    def set_halo(self, cells):
        """
        Update the halo
        :param cells:
        :return:
        """
        self.halo_cells = []
        for i in range(len(cells)):
            self.halo_cells.append(cells[i])

    def get(self, address, location='unknown'):
        """

        :param address: Address of cell to return
        :param location: Optional location value ("grid" or "halo")
        :return:
        """

        if location == 'unknown':
            if self.address_is_on_grid(address):
                location = 'grid'
            else:
                location = 'halo'

        if location == 'grid':
            address = tuple(address)
            return self.grid[address]
        elif location == 'halo':
            try:
                index = self.halo_addresses.index(address)
                return self.halo_cells[index]
            except ValueError:
                raise Exception("Address {0} is not on grid or in halo".format(address))

    def get_attribute(self, address, attribute, location='unknown'):
        """
        Get an attribute from a cell on the active grid (or halo)
        :param address:
        :param attribute:
        :param location: Optional location ("grid" or "halo")
        :return:
        """
        cell = self.get(address, location)
        if cell is None:
            return None
        else:
            return cell[attribute]

    def set_attribute_grid(self, address, attribute, value):
        """
        Set an attribute of a cell on the active grid (for initialisation and pre-processing)
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
        if attribute in self.attributes:
            self.work_grid[address][attribute] = value
        else:  # Specified attribute hasn't been set as a possibility
            raise Exception('Attribute {0} does not exist'.format(attribute))


class Neighbourhood:

    def __init__(self, dimensions, max_depth):
        self.dimensions = dimensions
        self.neighbour_table_moore, self.neighbour_table_von_neumann = self.construct_neighbour_table(max_depth)

    def construct_neighbour_table(self, max_depth):
        """
        Creates the required tables listing the addresses for varying depths for Moore and Von Neumann neighbourhoods.
        This is done once at the start so that we don't keep dynamically generating addresses each time.
        :param max_depth:
        :return:
        """
        moore_table = dict()
        von_neumann_table = dict()

        # Set up the required entries for Von Neumann
        for d in range(1, max_depth + 1):
            von_neumann_table[d] = []

        for d in range(max_depth):
            depth = d+1
            # Get truth table values
            range_ = range(-depth, depth + 1)
            row = list(itertools.product(range_, repeat=self.dimensions))
            # Remove the (0,0) entry
            row.remove((0,) * self.dimensions)
            reduced_row_moore = []
            von_neumann_table[depth] = []
            for neighbour in row:
                # Calculate Manhattan distance and add to appropraite von Neumann table row
                manhattan_distance = int(sum([math.fabs(x) for x in neighbour]))
                if manhattan_distance <= max_depth and neighbour not in von_neumann_table[manhattan_distance]:
                    von_neumann_table[manhattan_distance].append(neighbour)
                # Check if one of coordinates = depth, if so then use for moore at this depth
                for x in neighbour:
                    if int(math.fabs(x)) == depth:
                        reduced_row_moore.append(neighbour)
                        break
            moore_table[depth] = reduced_row_moore
        return moore_table, von_neumann_table

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
            output.append([address[j] + table[i][j] for j in range(len(address))])
        return output

    def neighbours_moore(self, address, depth=1):
        """
        Address of neighbours in the Moore neighbourhood of given depth
        :param address:
        :param depth:
        :param inclusive: Include lower depths?
        :return:
        """
        table = self.neighbour_table_moore[depth]
        return self.calculate_neighbours_locations(address, table)

    def neighbours_von_neumann(self, address, depth=1):
        """
        Gives the neighbours in the Von Neumann neighbourhood of address of given depth
        :param address:
        :param depth:
        :param inclusive: Includes lower depths?
        :return:
        """
        table = self.neighbour_table_von_neumann[depth]
        return self.calculate_neighbours_locations(address, table)


class EventHandler:

    def __init__(self):
        pass

    def handle_event(self, event):
        """
        Given an event, call the appropriate function to make changes to the grid
        :param event:
        :return:
        """
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
            raise Exception("Event ", type(event), "not handled")

    def process_bacterium_replication(self, event):
        """
        Add a new bacteria to the grid
        :param event:
        :return:
        """
        print "BACTERIUM REPLICATION"
        # Only process if the new bacterium address is on the grid
        if self.address_is_on_grid(event.new_bacterium_address):
            self.add_bacterium(event.new_bacterium_address, event.new_metabolism)

        if self.address_is_on_grid(event.original_bacterium_address):
            bacterium = self.get_attribute(event.original_bacterium_address, 'contents')
            bacterium.age = 0.0
            # Swap the division neighbourhood
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
        print "T CELL RECRUITMENT"
        # Only process if address is on the grid
        if self.address_is_on_grid(event.t_cell_address):
            self.add_t_cell(event.t_cell_address)

    def process_macrophage_recruitment(self, event):
        """
        Recruit a new macrophage from blood_vessel
        :param event:
        :return:
        """
        print "MACROPHAGE RECRUITMENT"
        if self.address_is_on_grid(event.macrophage_address):
            self.add_macrophage(event.macrophage_address, "resting")

    def process_chemo_kill_bacterium(self, event):
        """
        Chemotherapy kills a bacterium
        :param event:
        :return:
        """
        print "CHEMO KILL BACTERIUM"
        bacterium = self.get_attribute(event.address, 'contents')
        self.bacteria.remove(bacterium)
        #self.set_attribute_work_grid(event.address, 'contents', 0.0)

    def process_chemo_kill_macrophage(self, event):
        """
        Chemotherapy kills an infected macrophage
        :param event:
        :return:
        """
        print "CHEMO KILL MACROPHAGE"
        macrophage = self.get_attribute(event.dependant_addresses[0], 'contents')
        self.macrophages.remove(macrophage)
        #self.set_attribute_work_grid(event.macrophage_to_kill.address, 'contents', 0.0)

    def process_t_cell_death(self, event):
        """
        A t-cell dies (through age)
        :param event:
        :return:
        """
        print "T-CELL DEATH"
        t_cell_to_die = self.get_attribute(event.address, 'contents')
        #self.set_attribute_work_grid(event.address, 'contents', 0.0)
        self.t_cells.remove(t_cell_to_die)

    def process_t_cell_movement(self, event):
        """
        A t-cell moves to a new cell
        :param event:
        :return:
        """
        print "T-CELL MOVEMENT"
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        # T-cell moving between 2 cells in the same tile
        if event.internal:
            t_cell = self.get_attribute(from_address, 'contents')
            t_cell.address = to_address
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            #self.set_attribute_work_grid(to_address, 'contents', t_cell)
        elif self.address_is_on_grid(from_address):  # T-cell is moving to a new tile
            t_cell = self.get_attribute(from_address, 'contents')
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            self.t_cells.remove(t_cell)
        elif self.address_is_on_grid(to_address):  # T-cell has arrived from another tile
            event.t_cell_to_move.address = to_address
            #self.set_attribute_work_grid(to_address, 'contents', event.t_cell_to_move)
            self.t_cells.append(event.t_cell_to_move)

    def process_t_cell_kill_macrophage(self, event):
        print "T-CELL KILLS MACROPHAGE"
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        if self.address_is_on_grid(to_address):
            # Turn macrophage into caseum
            macrophage = self.get_attribute(to_address, 'contents')
            self.macrophages.remove(macrophage)
            # self.set_attribute_work_grid(to_address, 'contents', 'caseum')
            self.caseum.append(to_address)

        if self.address_is_on_grid(from_address):
            # Remove t-cell
            t_cell = self.get_attribute(from_address, 'contents')
            self.t_cells.remove(t_cell)
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)

    def process_macrophage_death(self, event):
        print "MACROPHAGE DEATH"
        # Resting or active die, infected/chronically infected turn to caseum
        macrophage_to_die = self.get_attribute(event.address, 'contents')
        #if macrophage_to_die.state == 'resting' or macrophage_to_die.state == 'active':
        #    self.set_attribute_work_grid(macrophage_to_die.address, 'contents', 0.0)
        #el
        if macrophage_to_die.state == 'infected' or macrophage_to_die.state == 'chronically_infected':
            # self.set_attribute_work_grid(macrophage_to_die.address, 'contents', 'caseum')
            self.caseum.append(macrophage_to_die.address)
        # Remove macrophage
        self.macrophages.remove(macrophage_to_die)

    def process_macrophage_movement(self, event):
        print "MACROPHAGE MOVEMENT"
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]
        # Macrophage moving between 2 cells in the same tile
        if event.internal:
            macrophage = self.get_attribute(from_address, 'contents')
            macrophage.address = to_address
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            #self.set_attribute_work_grid(to_address, 'contents', macrophage)
        elif self.address_is_on_grid(from_address):  # Macrophage is moving to a new tile
            # Remove macrophage
            macrophage = self.get_attribute(from_address, 'contents')
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            self.macrophages.remove(macrophage)
        elif self.address_is_on_grid(to_address):  # Macrophage has arrived from another tile
            # Add macrophage
            event.macrophage_to_move.address = to_address
            #self.set_attribute_work_grid(to_address, 'contents', event.macrophage_to_move)
            self.macrophages.append(event.macrophage_to_move)

    def process_macrophage_kills_bacterium(self, event):
        print "MACROPHAGE_KILLS_BACTERIUM"
        from_address = event.dependant_addresses[0]
        to_address = event.dependant_addresses[1]

        # Macrophage moving between 2 cells in the same tile
        if event.internal:
            # Move macrophage
            macrophage = self.get_attribute(from_address, 'contents')
            macrophage.address = to_address
            # Remove bacterium
            bacterium = self.get_attribute(to_address, 'contents')
            self.bacteria.remove(bacterium)
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            #self.set_attribute_work_grid(to_address, 'contents', macrophage)

            if macrophage.state == 'resting' or macrophage.state == 'infected' or macrophage.state == \
                    'chronically_infected':
                # Macrophage ingests bacteria, doesn't kill
                event.macrophage_to_move.intracellular_bacteria += 1
        elif self.address_is_on_grid(from_address):  # Macrophage is moving to a new tile
            macrophage = self.get_attribute(from_address, 'contents')
            #self.set_attribute_work_grid(from_address, 'contents', 0.0)
            self.macrophages.remove(macrophage)
        elif self.address_is_on_grid(to_address):  # Macrophage has arrived from another tile
            event.macrophage_to_move.address = to_address
            bacterium = self.get_attribute(to_address, 'contents')
            self.bacteria.remove(bacterium)
            #self.set_attribute_work_grid(to_address, 'contents', event.macrophage_to_move)
            self.macrophages.append(event.macrophage_to_move)

            if event.macrophage_to_move.state == 'resting' or event.macrophage_to_move.state == 'infected' or \
                    event.macrophage_to_move.state == 'chronically_infected':
                # Macrophage ingests bacterium, doesn't kill
                event.macrophage_to_move.intracellular_bacteria += 1

    def process_macrophage_state_change(self, event):
        print "MACROPHAGE_STATE_CHANGE: to", event.new_state
        # Pulling the macrophage from the grid
        macrophage = self.get_attribute(event.address, 'contents')
        macrophage.state = event.new_state
        #self.set_attribute_work_grid(event.address, 'contents', macrophage)

    def process_bacterium_state_change(self, event):
        print "BACTERIUM_STATE_CHANGE:", event.type_of_change, " to", event.new_value
        # Pull bacterium from grid
        bacterium = self.get_attribute(event.address, 'contents')

        if event.type_of_change == 'metabolism':
            bacterium.metabolism = event.new_value
        elif event.type_of_change == 'resting':
            bacterium.resting = event.new_value
        #self.set_attribute_work_grid(event.address, 'contents', bacterium)

    def process_macrophage_bursting(self, event):
        print "MACROPHAGE BURSTING"
        macrophage_to_burst = self.get_attribute(event.macrophage_address, 'contents')
        self.set_attribute_work_grid(macrophage_to_burst.address, 'contents', 'caseum')
        self.caseum.append(macrophage_to_burst.address)
        self.macrophages.remove(macrophage_to_burst)

        # TODO - COMP - bacteria stuff
        for i in event.bacteria_addresses:
            if i in event.impacted_addresses_allowed and self.address_is_on_grid(i):
                    self.add_bacterium(i, 'slow')


class Automaton(Tile, Neighbourhood, EventHandler):

    def __init__(self, shape, tile_id, attributes, parameters, blood_vessels, fast_bacteria=None, slow_bacteria=None,
                 macrophages=None):
        Tile.__init__(self, shape, attributes)
        Neighbourhood.__init__(self, len(shape), parameters['max_depth'])
        EventHandler.__init__(self)

        # Max depth must be +1 or more greater than the caseum distance as we need to work out diffusion rates one
        # cell deep into the halo
        assert parameters['caseum_distance_to_reduce_diffusion'] + 1 <= parameters['max_depth']

        self.tile_id = tile_id
        self.parameters = parameters
        self.time = 0

        self.max_oxygen_local = 0.0
        self.max_chemotherapy_local = 0.0
        self.max_chemokine_local = 0.0

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

        # INITIAL VESSELS
        self.initialise_blood_vessels(blood_vessels)
        self.initialise_oxygen_levels()

        # COPY GRID TO WORK GRID
        self.create_work_grid()

        # INITIAL AGENTS
        self.initialise_bacteria(fast_bacteria, slow_bacteria)
        self.initialise_macrophages(macrophages)

        self.swap_grids()

    def initialise_blood_vessels(self, addresses):
        for address in addresses:
            self.set_attribute_grid(address, 'blood_vessel', 1.5)
            self.blood_vessels.append(address)

    def initialise_bacteria(self, fast_bacteria, slow_bacteria):
        for address in fast_bacteria:
            self.add_bacterium(address, "fast")
        for address in slow_bacteria:
            self.add_bacterium(address, "slow")

    def initialise_macrophages(self, addresses):
        for address in addresses:
            self.add_macrophage(address, "resting")

    def initialise_oxygen_levels(self):
        for address in self.blood_vessels:
            self.set_attribute_grid(address, 'oxygen', self.parameters['initial_oxygen'])
            self.max_oxygen_local = max(self.max_oxygen_local, self.parameters['initial_oxygen'])

    def update(self):
        """
        Performs the continuous update (diffusion of oxygen, chemokine, chemotherapy) and applies to the work grid,
        and creates a list of potential events to be executed.
        Events are based on the state of the current GRID (i.e. are not affected by the diffusion changes)
        :return:
        """

        self.time += 1

        # ----------------------------
        # CONTINUOUS (Diffusion)
        # ----------------------------

        # Pre-processing (calculating diffusion rates)
        self.diffusion_pre_process_v2()
        # In chemo window?
        if (((self.parameters['chemotherapy_schedule1_start'] / self.parameters['time_step']) <=
                self.time <
                (self.parameters['chemotherapy_schedule1_end'] / self.parameters['time_step'])) or
                (self.parameters['chemotherapy_schedule2_start'] / self.parameters['time_step'] <= self.time)):
            chemo = True
        else:
            chemo = False

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

        # BACTERIA REPLICATION
        # If bacteria is of suitable age, attempts to replicate (create a duplicate in neighbouring cell)
        self.bacteria_replication()

        # T-CELL RECRUITMENT
        self.t_cell_recruitment()

        # MACROPHAGE RECRUITMENT
        self.macrophage_recruitment()

        # CHEMOTHERAPY KILLING BACTERIA
        self.chemotherapy_killing_bacteria()

        # CHEMOTHERAPY KILLING MACROPHAGES
        self.chemotherapy_killing_macrophages()

        # T-CELL DEATH, MOVEMENT & MACROPHAGE KILLING
        self.t_cell_processes()

        # MACROPHAGES - death, movement, bacteria ingestion
        self.macrophage_processes()

        # MACROPHAGE STATE CHANGES / BURSTING
        self.macrophage_state_changes()

        # BACTERIUM STATE CHANGES
        self.bacteria_state_changes()

        # Reorder events
        self.reorder_events()

    def reorder_events(self):
        # TODO - COMP - other methods - currently just random
        np.random.shuffle(self.potential_events)

    def process_events(self, events):
        """
        Given a list of acceptable addresses, pass to the Event Handler to update the grid
        :param events:
        :return:
        """
        for event in events:
            self.handle_event(event)

        # Ensure all agents are put onto the work grid
        self.persist_agents()
        # Swap work grid with main grid
        self.swap_grids()

    def diffusion_pre_process(self):
        """
        Pre-processing for diffusion. Looks through all cells and calculates the diffusion rates - diffusion rate
        of oxygen and chemotherapy drops if there is too much caseum in the vicinity
        :return:
        """
        # Loop through all cells
        for location in range(self.size):
            address = self.location_to_address(location)

            # Get initial diffusion rates
            oxygen_diffusion = self.parameters['oxygen_diffusion']
            chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

            # Get all neighbours up to the specified distance
            neighbours = []
            for i in range(1, int(self.parameters['caseum_distance_to_reduce_diffusion']) + 1):
                neighbours += self.neighbours_moore(address, i)

            # Check if there is specified amount of caseum within specified distance of cell
            caseum_count = 0
            for neighbour_address in neighbours:
                cell = self.get(neighbour_address)
                if cell is not None and cell['contents'] == 'caseum':
                    caseum_count += 1
                    # Once the caseum threshold is reached
                    if caseum_count == self.parameters['caseum_threshold_to_reduce_diffusion']:
                        # Decrease the diffusion level at the cell
                        oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                        chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                        # Exit the loop
                        break

            # Need to set the values on the current grid
            self.set_attribute_grid(address, 'oxygen_diffusion_rate', oxygen_diffusion)
            self.set_attribute_grid(address, 'chemotherapy_diffusion_rate', chemotherapy_diffusion)

        # Set diffusion rates on halo depth 1
        for halo_address in self.halo_depth1:
            if self.get(halo_address, "halo") is not None:

                # Get initial rates
                oxygen_diffusion = self.parameters['oxygen_diffusion']
                chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']
                # Get neighbours
                neighbours = []
                for i in range(1, int(self.parameters['caseum_distance_to_reduce_diffusion']) + 1):
                    neighbours += self.neighbours_moore(halo_address, i)
                # Check if caseum in neighbourhood exceeds threshold
                caseum_count = 0
                for neighbour_address in neighbours:
                    cell = self.get(neighbour_address)
                    if cell is not None and cell['contents'] == 'caseum':
                        caseum_count += 1
                        # Once the caseum threshold is reached
                        if caseum_count == self.parameters['caseum_threshold_to_reduce_diffusion']:
                            # Decrease the diffusion level at the cell
                            oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                            chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                            # Exit the loop
                            break

                # Need to set the values on the halo
                index = self.halo_addresses.index(halo_address)
                self.halo_cells[index]['oxygen_diffusion_rate'] = oxygen_diffusion
                self.halo_cells[index]['chemotherapy_diffusion_rate'] = chemotherapy_diffusion

    def diffusion_pre_process_v2(self):

        affected_addresses = []
        caseum_addresses = list(self.caseum)
        
        for index in range(len(self.halo_cells)):
            if self.halo_cells[index] is not None and self.halo_cells[index]['contents'] == 'caseum':
                caseum_addresses.append(self.halo_addresses[index])

        for address in caseum_addresses:

            for depth in range(1, self.parameters['caseum_distance_to_reduce_diffusion']+1):
                neighbours = self.neighbours_moore(address, depth)
                for n in neighbours:
                    tuple_address = tuple(n)
                    affected_addresses.append(tuple_address)

        counted = Counter(affected_addresses)

        for location in range(self.size):
            address = self.location_to_address(location)

            # Get initial diffusion rates
            oxygen_diffusion = self.parameters['oxygen_diffusion']
            chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

            if tuple(address) in counted and counted[tuple(address)] >= \
                    self.parameters['caseum_threshold_to_reduce_diffusion']:
                oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']

            # Need to set the values on the current grid
            self.set_attribute_grid(address, 'oxygen_diffusion_rate', oxygen_diffusion)
            self.set_attribute_grid(address, 'chemotherapy_diffusion_rate', chemotherapy_diffusion)

        for halo_address in self.halo_depth1:
            if self.get(halo_address, "halo") is not None:

                # Get initial rates
                oxygen_diffusion = self.parameters['oxygen_diffusion']
                chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

                if tuple(halo_address) in counted and counted[tuple(halo_address)] >= \
                        self.parameters['caseum_threshold_to_reduce_diffusion']:
                    oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                    chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']

                # Need to set the values on the halo
                index = self.halo_addresses.index(halo_address)

                self.halo_cells[index]['oxygen_diffusion_rate'] = oxygen_diffusion
                self.halo_cells[index]['chemotherapy_diffusion_rate'] = chemotherapy_diffusion

    def diffusion(self, chemo):
        for i in range(self.size):
            address = self.location_to_address(i)

            # OXYGEN
            oxygen_level = self.oxygen(address)
            self.set_attribute_work_grid(address, 'oxygen', oxygen_level)

            # CHEMOTHERAPY
            if chemo:
                chemotherapy_level = self.chemotherapy(address)
                self.set_attribute_work_grid(address, 'chemotherapy', chemotherapy_level)
            else:
                # TODO - MED - check validity of this (how can chemotherapy suddenly disappear)
                self.set_attribute_work_grid(address, 'chemotherapy', 0.0)

            # CHEMOKINE
            chemokine_level = self.chemokine(address)
            self.set_attribute_work_grid(address, 'chemokine', chemokine_level)

    def bacteria_replication(self):
        for bacterium in self.bacteria:
            bacterium.age += self.parameters['time_step']

            # Skip if the bacterium is resting
            if bacterium.resting:
                continue

            division = False
            if bacterium.metabolism == 'fast':
                maximum = self.parameters['bacteria_replication_fast_upper']
                minimum = self.parameters['bacteria_replication_fast_lower']
            else:  # Slow
                maximum = self.parameters['bacteria_replication_slow_upper']
                minimum = self.parameters['bacteria_replication_slow_lower']

            replication_time = np.random.randint(minimum, maximum)

            if self.time % replication_time == 0:
                division = True

            if division:
                free_neighbours = []

                for depth in range(1, 4):
                    if bacterium.division_neighbourhood == 'mo':
                        neighbours = self.neighbours_moore(bacterium.address, depth)
                    else:
                        neighbours = self.neighbours_von_neumann(bacterium.address, depth)

                    for neighbour_address in neighbours:
                        neighbour = self.get(neighbour_address)
                        if neighbour is not None and neighbour['contents'] == 0.0 and \
                                        neighbour_address not in self.blood_vessels:
                            free_neighbours.append(neighbour_address)

                    if len(free_neighbours) > 0:
                        break

                if len(free_neighbours) == 0:
                    new_event = BacteriumStateChange(bacterium.address, 'resting', True)
                    self.potential_events.append(new_event)
                else:  # Free space found
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]

                    if self.address_is_on_grid(neighbour_address):
                        internal = True
                    else:
                        internal = False
                    new_event = BacteriumReplication(bacterium.address, neighbour_address, bacterium.metabolism,
                                                     internal)
                    self.potential_events.append(new_event)

    def t_cell_recruitment(self):
        if self.number_of_bacteria_global >= self.parameters['bacteria_threshold_for_t_cells']:
            # Each blood vessel
            for bv_address in self.blood_vessels:
                r = np.random.randint(1, 100)
                if r <= self.parameters['t_cell_recruitment_probability']:
                    neighbours = self.neighbours_von_neumann(bv_address, 1)
                    # Get neighbours which are empty and have a high enough chemokine level
                    free_neighbours = []
                    for neighbour_address in neighbours:
                        if self.get(neighbour_address) is not None and \
                                        self.get_attribute(neighbour_address, 'blood_vessel') == 0.0 and \
                                        self.get_attribute(neighbour_address, 'contents') == 0.0 and \
                                        self.chemokine_scale(neighbour_address) > self.parameters[
                                    'chemokine_scale_for_t_cell_recruitment']:
                            free_neighbours.append(neighbour_address)
                    # Check there is free space
                    if len(free_neighbours) > 0:
                        # Pick one of the neighbours
                        neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                        if self.address_is_on_grid(neighbour_address):
                            internal = True
                        else:
                            internal = False
                        new_event = RecruitTCell(neighbour_address, internal)
                        self.potential_events.append(new_event)

    def macrophage_recruitment(self):
        for bv_address in self.blood_vessels:
            r = np.random.randint(1, 100)
            if r <= self.parameters['macrophage_recruitment_probability']:
                neighbours = self.neighbours_von_neumann(bv_address, 1)
                free_neighbours = []
                for neighbour_address in neighbours:
                    if self.get(neighbour_address) is not None and \
                                    self.get_attribute(neighbour_address, 'blood_vessel') == 0.0 and \
                                    self.get_attribute(neighbour_address, 'contents') == 0.0 and \
                                    self.chemokine_scale(neighbour_address) > self.parameters[
                                'chemokine_scale_for_macrophage_recruitment']:
                        free_neighbours.append(neighbour_address)

                if len(free_neighbours) > 0:
                    # Pick one of the neighbours
                    chosen_neighbour = free_neighbours[np.random.randint(len(free_neighbours))]
                    if self.address_is_on_grid(chosen_neighbour):
                        internal = True
                    else:
                        internal = False
                    new_event = RecruitMacrophage(chosen_neighbour, internal)
                    self.potential_events.append(new_event)

    def chemotherapy_killing_bacteria(self):
        for bacterium in self.bacteria:
            chemo_scale = self.chemotherapy_scale(bacterium.address)
            if (bacterium.metabolism == 'fast' and chemo_scale >
                self.parameters['chemotherapy_scale_for_kill_fast_bacteria']) \
                    or \
                    (bacterium.metabolism == 'slow' and chemo_scale >
                        self.parameters['chemotherapy_scale_for_kill_slow_bacteria']):
                new_event = ChemoKillBacterium(bacterium.address)
                self.potential_events.append(new_event)

    def chemotherapy_killing_macrophages(self):
        for m in self.macrophages:
            chemo_scale = self.chemotherapy_scale(m.address)
            if ((m.state == 'infected' or m.state == 'chronically_infected') and chemo_scale >
                self.parameters['chemotherapy_scale_for_kill_macrophage']):
                new_event = ChemoKillMacrophage(m.address)
                self.potential_events.append(new_event)

    def t_cell_processes(self):
        # TODO - MED - does this make sense - e.g. t-cell death is dependent on the bacteria number and movement time
        if self.number_of_bacteria_global >= self.parameters['bacteria_threshold_for_t_cells'] and \
                                self.time % self.parameters['t_cell_movement_time'] == 0:

            for t_cell in self.t_cells:
                t_cell.age += self.parameters['time_step']

                age_threshold = np.random.randint(0, self.parameters['t_cell_age_threshold'])

                # T-CELL DEATH
                if t_cell.age >= age_threshold:
                    new_event = TCellDeath(t_cell.address)
                    self.potential_events.append(new_event)
                else:  # T-CELL MOVE
                    random_move = False
                    prob_random_move = np.random.randint(1, 101)
                    if prob_random_move <= self.parameters['t_cell_random_move_probability']:
                        random_move = True

                    neighbours = self.neighbours_moore(t_cell.address, 1)
                    if random_move:
                        possible_neighbours = []
                        for n in neighbours:
                            if self.get(n) is not None:
                                possible_neighbours.append(n)
                        index = np.random.randint(0, len(possible_neighbours))
                        chosen_neighbour_address = possible_neighbours[index]
                    else:
                        chosen_index = self.find_max_chemokine_neighbour(neighbours)[0]
                        chosen_neighbour_address = neighbours[chosen_index]

                    # Check if leaving the grid
                    internal = self.address_is_on_grid(chosen_neighbour_address)

                    if internal:
                        location = 'grid'
                    else:
                        location = 'halo'

                    neighbour = self.get(chosen_neighbour_address, location)

                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = TCellMovement(t_cell, t_cell.address, chosen_neighbour_address, internal)
                        self.potential_events.append(new_event)
                    elif isinstance(neighbour['contents'], Macrophage) and (neighbour['contents'].state == 'infected'
                                                                            or neighbour[
                            'contents'].state == 'chronically_infected'):
                        prob_tcell_killing = np.random.randint(1, 101)
                        if prob_tcell_killing <= self.parameters['t_cell_kills_macrophage_probability']:
                            new_event = TCellKillsMacrophage(t_cell, t_cell.address, chosen_neighbour_address, internal)
                            self.potential_events.append(new_event)

    def macrophage_processes(self):
        for macrophage in self.macrophages:
            macrophage.age += self.parameters['time_step']

            if macrophage.state == 'resting':

                # TODO - COMP/MED - can't do 0 as get division errors
                random_macrophage_age = np.random.randint(1, self.parameters['resting_macrophage_age_limit'])

                if macrophage.age % random_macrophage_age == 0:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue

                # TODO - MED - movement based on time not age
                if self.time % self.parameters['resting_macrophage_movement_time'] == 0:

                    neighbours = self.neighbours_moore(macrophage.address, 1)
                    chosen_index, max_chemokine_scale = self.find_max_chemokine_neighbour(neighbours)

                    prob_random_move = np.random.randint(1, 101)
                    random_move = False
                    if prob_random_move <= self.parameters['prob_resting_macrophage_random_move'] \
                            or max_chemokine_scale <= \
                                    self.parameters['minimum_chemokine_for_resting_macrophage_movement']:
                        random_move = True

                    if random_move:
                        chosen_neighbour_address = neighbours[np.random.randint(0, len(neighbours))]
                    else:
                        chosen_neighbour_address = neighbours[chosen_index]

                    # Check if leaving the grid
                    internal = self.address_is_on_grid(chosen_neighbour_address)

                    if internal:
                        location = 'grid'
                    else:
                        location = 'halo'

                    neighbour = self.get(chosen_neighbour_address,location)

                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)

            elif macrophage.state == 'active':
                # TODO - MED - death is based on age > limit, no prob
                if macrophage.age > self.parameters['active_macrophage_age_limit']:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue

                if self.time % self.parameters['active_macrophage_movement_time'] == 0:
                    neighbours = self.neighbours_moore(macrophage.address, 1)
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]

                    internal = self.address_is_on_grid(chosen_neighbour_address)
                    if internal:
                        location = 'grid'
                    else:
                        location = 'halo'
                    neighbour = self.get(chosen_neighbour_address, location)

                    if isinstance(neighbour['contents'], Bacterium):

                        prob_macrophage_kill = np.random.randint(1, 101)
                        if (neighbour['contents'].metabolism == 'fast' and prob_macrophage_kill <= self.parameters[
                            'prob_active_macrophage_kill_fast_bacteria']) or (
                                        neighbour['contents'].metabolism == 'slow' and prob_macrophage_kill <=
                                    self.parameters[
                                        'prob_active_macrophage_kill_slow_bacteria']):
                            new_event = MacrophageKillsBacterium(macrophage, macrophage.address,
                                                                 chosen_neighbour_address, internal)
                            self.potential_events.append(new_event)

                    elif neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)

            elif macrophage.state == 'infected':
                # TODO - COMP/MED - can't do 0 as get division errors
                random_macrophage_age = np.random.randint(1, self.parameters['infected_macrophage_age_limit'])

                if macrophage.age % random_macrophage_age == 0:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue

                if self.time % self.parameters['infected_macrophage_movement_time'] == 0:
                    neighbours = self.neighbours_moore(macrophage.address, 1)
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]

                    internal = False
                    location = 'halo'
                    if self.address_is_on_grid(chosen_neighbour_address):
                        internal = True
                        location = 'grid'

                    neighbour = self.get(chosen_neighbour_address, location)

                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)

            elif macrophage.state == 'chronically_infected':
                # TODO - COMP/MED - can't do 0 as get division errors

                random_macrophage_age = np.random.randint(1,
                                                          self.parameters['chronically_infected_macrophage_age_limit'])

                if macrophage.age % random_macrophage_age == 0:
                    new_event = MacrophageDeath(macrophage.address)
                    self.potential_events.append(new_event)
                    # Progress to the next macrophage
                    continue

                if self.time % self.parameters['chronically_infected_macrophage_movement_time'] == 0:
                    neighbours = self.neighbours_moore(macrophage.address, 1)
                    chosen_neighbour_address = neighbours[self.find_max_chemokine_neighbour(neighbours)[0]]

                    internal = False
                    location = 'halo'
                    if self.address_is_on_grid(chosen_neighbour_address):
                        internal = True
                        location = 'grid'

                    neighbour = self.get(chosen_neighbour_address, location)

                    if neighbour['contents'] == 0.0 and neighbour['blood_vessel'] == 0.0:
                        new_event = MacrophageMovement(macrophage, macrophage.address, chosen_neighbour_address,
                                                       internal)
                        self.potential_events.append(new_event)
                    elif isinstance(neighbour['contents'], Bacterium):
                        new_event = MacrophageKillsBacterium(macrophage, macrophage.address, chosen_neighbour_address,
                                                             internal)
                        self.potential_events.append(new_event)

    def macrophage_state_changes(self):
        for macrophage in self.macrophages:

            # MACROPHAGE STATE CHANGES
            if macrophage.state == 'resting':
                # Resting to active
                if self.chemokine_scale(macrophage.address) > \
                        self.parameters['chemokine_scale_for_macrophage_activation'] and \
                                macrophage.intracellular_bacteria == 0:
                    new_event = MacrophageChangesState(macrophage.address, "active")
                    self.potential_events.append(new_event)
                # Resting to infected
                elif macrophage.intracellular_bacteria == 1:
                    new_event = MacrophageChangesState(macrophage.address, "infected")
                    self.potential_events.append(new_event)
            elif macrophage.state == 'active':
                # Active to resting
                if self.chemokine_scale(macrophage.address) < \
                        self.parameters['chemokine_scale_for_macrophage_deactivation']:
                    new_event = MacrophageChangesState(macrophage.address, "resting")
                    self.potential_events.append(new_event)
            elif macrophage.state == 'infected':
                # Infected to Chronically Infected
                if macrophage.intracellular_bacteria > self.parameters['bacteria_to_turn_chronically_infected']:
                    new_event = MacrophageChangesState(macrophage.address, "chronically_infected")
                    self.potential_events.append(new_event)
            elif macrophage.state == 'chronically_infected':
                # MACROPHAGE BURSTING
                if macrophage.intracellular_bacteria == self.parameters['bacteria_to_burst_macrophage']:
                    internal = self.address_is_on_grid(macrophage.address)
                    bacteria_addresses = []
                    for depth in range(1, 4):
                        neighbours = self.neighbours_moore(macrophage.address, depth)
                        # Shuffle the neighbours so we don't give priority
                        np.random.shuffle(neighbours)
                        for n in neighbours:
                            # TODO - COMP - do this (contents and BV) too often.
                            # Could maybe make blood vessel part of contents
                            if self.get(n) is not None and self.get_attribute(n, 'contents') == 0.0 and \
                                            self.get_attribute(n, 'blood_vessel') == 0.0:
                                bacteria_addresses.append(n)
                                if not self.address_is_on_grid(n):
                                    internal = False
                            if len(bacteria_addresses) == self.parameters['bacteria_to_burst_macrophage']:
                                break
                        if len(bacteria_addresses) == self.parameters['bacteria_to_burst_macrophage']:
                            break
                    new_event = MacrophageBursting(macrophage.address, bacteria_addresses, internal)
                    self.potential_events.append(new_event)

    def bacteria_state_changes(self):
        for bacterium in self.bacteria:

            # Metabolism change only happens later in process
            if self.time > 2 / self.parameters['time_step']:

                if bacterium.metabolism == 'fast' and self.oxygen_scale(bacterium.address) <= self.parameters[
                    'oxygen_scale_for_metabolism_change_to_slow']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'slow')
                    self.potential_events.append(new_event)
                elif bacterium.metabolism == 'slow' and self.oxygen_scale(bacterium.address) > self.parameters[
                    'oxygen_scale_for_metabolism_change_to_fast']:
                    new_event = BacteriumStateChange(bacterium.address, 'metabolism', 'fast')
                    self.potential_events.append(new_event)

            if bacterium.resting:
                space_found = False
                for depth in range(1, 4):
                    neighbours = self.neighbours_moore(bacterium.address, depth)
                    for n in neighbours:
                        if self.get(n) is not None and self.get_attribute(n, 'blood_vessel') == 0.0 and \
                                        self.get_attribute(n, 'contents') == 0.0:
                            new_event = BacteriumStateChange(bacterium.address, 'resting', False)
                            self.potential_events.append(new_event)
                            space_found = True
                            break
                    if space_found:
                        break

    def oxygen(self, address):
        # Get the current cell values
        cell = self.get(address, "grid")

        # Get diffusion value for cell
        cell_diffusion = cell['oxygen_diffusion_rate']

        # Initialise expression
        expression = 0

        # Get immediate von neumann neighbours
        neighbour_addresses = self.neighbours_von_neumann(address, 1)
        for neighbour_address in neighbour_addresses:
            neighbour = self.get(neighbour_address)
            # Only process if not boundary
            if neighbour is not None:
                neighbour_diffusion = neighbour['oxygen_diffusion_rate']
                expression += ((cell_diffusion + neighbour_diffusion) / 2 * (neighbour['oxygen'] - cell['oxygen'])) / \
                              (self.parameters['spatial_step'] * self.parameters['spatial_step'])

        # Add oxygen entering through blood vessel
        expression += (self.parameters['oxygen_from_source'] * cell['blood_vessel'])
        # If there is bacteria in cell, then oxygen is taken up by bacteria so remove
        if isinstance(cell['contents'], Bacterium):
            expression -= self.parameters['oxygen_uptake_from_bacteria'] * cell['oxygen']

        # Calculate new level
        new_oxygen = cell['oxygen'] + self.parameters['time_step'] * expression

        # Overwrite the maximum oxygen value if larger
        self.max_oxygen_local = max(self.max_oxygen_local, new_oxygen)

        return new_oxygen

    def chemotherapy(self, address):
        # Get the current cell values
        cell = self.get(address, "grid")

        # Get diffusion value for cell based on system parameters
        cell_diffusion = cell['chemotherapy_diffusion_rate']

        # Initialise expression
        expression = 0

        # Get immediate von Neumann neighbours
        neighbour_addresses = self.neighbours_von_neumann(address, 1)
        for neighbour_address in neighbour_addresses:
            neighbour = self.get(neighbour_address)
            # Only process if not boundary
            if neighbour is not None:
                # Pre-calculated diffusion rate
                neighbour_diffusion = neighbour['chemotherapy_diffusion_rate']

                expression += ((cell_diffusion + neighbour_diffusion) / 2 * (
                    neighbour['chemotherapy'] - cell['chemotherapy'])) / (
                                  self.parameters['spatial_step'] * self.parameters['spatial_step'])

        # Release of chemotherapy from blood vessel
        expression += (self.parameters['chemotherapy_from_source'] * cell['blood_vessel'])

        # Chemotherapy decay
        expression -= self.parameters['chemotherapy_decay'] * cell['chemotherapy']

        # Calculate new level
        new_chemotherapy = cell['chemotherapy'] + self.parameters['time_step'] * expression

        self.max_chemotherapy_local = max(self.max_chemotherapy_local, new_chemotherapy)

        return new_chemotherapy

    def chemokine(self, address):

        # Get the current cell values
        cell = self.get(address, "grid")

        # Get diffusion value for cell based on system parameters
        cell_diffusion = self.parameters['chemokine_diffusion']

        # Initialise expression
        expression = 0

        # Get immediate von Neumann neighbours
        neighbour_addresses = self.neighbours_von_neumann(address, 1)

        for neighbour_address in neighbour_addresses:
            neighbour = self.get(neighbour_address)
            # Only process if not boundary
            if neighbour is not None:
                neighbour_diffusion = self.parameters['chemokine_diffusion']
                expression += ((cell_diffusion + neighbour_diffusion) / 2 * (
                    neighbour['chemokine'] - cell['chemokine'])) / (
                                  self.parameters['spatial_step'] * self.parameters['spatial_step'])

        # Release of chemokine by bacteria
        if isinstance(cell['contents'], Bacterium):
            expression += self.parameters['chemokine_from_bacteria']

        # Release of chemokine by macrophages
        if isinstance(cell['contents'], Macrophage):
            expression += self.parameters['chemokine_from_macrophage']

        # Chemokine decay
        expression -= self.parameters['chemokine_decay'] * cell['chemokine']

        # Calculate new level
        new_chemokine = cell['chemokine'] + self.parameters['time_step'] * expression

        self.max_chemokine_local = max(self.max_chemokine_local, new_chemokine)
        return new_chemokine

    def set_max_oxygen_global(self, max_oxygen):
        self.max_oxygen_global = max_oxygen

    def set_max_chemotherapy_global(self, max_chemotherapy):
        self.max_chemotherapy_global = max_chemotherapy

    def set_max_chemokine_global(self, max_chemokine):
        self.max_chemokine_global = max_chemokine

    def set_global_bacteria_number(self, number):
        self.number_of_bacteria_global = number

    def oxygen_scale(self, address):
        if self.max_oxygen_global == 0.0:
            return 0.0
        else:
            return (self.get_attribute(address, 'oxygen') / self.max_oxygen_global) * 100

    def chemotherapy_scale(self, address):
        if self.max_chemotherapy_global == 0.0:
            return 0.0
        else:
            return (self.get_attribute(address, 'chemotherapy') / self.max_chemotherapy_global) * 100

    def chemokine_scale(self, address):
        if self.max_chemokine_global == 0.0:
            return 0.0
        else:
            return (self.get_attribute(address, 'chemokine') / self.max_chemokine_global) * 100

    def add_bacterium(self, address, metabolism):
        new_bacterium = Bacterium(address, metabolism)
        self.bacteria.append(new_bacterium)
        self.set_attribute_work_grid(address, 'contents', new_bacterium)

    def add_macrophage(self, address, state):
        new_macrophage = Macrophage(address, state)
        self.macrophages.append(new_macrophage)
        self.set_attribute_work_grid(address, 'contents', new_macrophage)

    def add_t_cell(self, address):
        new_t_cell = TCell(address)
        self.t_cells.append(new_t_cell)
        self.set_attribute_work_grid(address, 'contents', new_t_cell)

    def persist_agents(self):
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
            if self.get(neighbours[index]) is not None:
                chemokine_scale = self.chemokine_scale(neighbours[index])
                if chemokine_scale > max_chemokine_scale:
                    max_chemokine_scale = chemokine_scale
                    highest_indices = [index]
                elif chemokine_scale == max_chemokine_scale:
                    highest_indices.append(index)

        # Tie-breaking. If just one pick it, else pick any one index at random
        if len(highest_indices) == 1:
            chosen_index = highest_indices[0]
        else:
            choice = np.random.randint(0, len(highest_indices))
            chosen_index = highest_indices[choice]

        return chosen_index, max_chemokine_scale

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