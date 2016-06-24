import numpy as np
import math
import itertools


class Topology:

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_local=[],
                 fast_bacteria_local=[], slow_bacteria_local=[], macrophages_local = []):

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
        # TODO - is there a better way of doing this?
        depth1_addresses = self.get_external_addresses_required(self.automata[0],1)

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
            neighbours = automaton.neighbours_moore(address, depth)
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

    def __init__(self, tile_arrangement, total_shape, attributes, parameters, blood_vessel_global=[],
                 fast_bacteria_global=[], slow_bacteria_global=[], macrophages_global = []):
        # Check it's two dimensions
        assert len(total_shape) == 2
        self.number_of_tiles = reduce(lambda x, y: x * y, tile_arrangement)
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

    def address_to_location(self, address):
        """
        Convert an address set of coordinates into an integer location
        :param address:
        :return:
        """
        # TODO - probably redundant function now
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
        """
        Get the cell values of cells in the danger zone
        :return:
        """
        danger_zone = []
        for address in self.danger_zone_addresses:
            # Get cell from the work grid
            danger_zone.append(self.work_grid[tuple(address)])

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

    def get(self, address):
        """
        Get a cell from the active grid (or from halo)
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
        Get an attribute from a cell on the active grid (or halo)
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
        # TODO - probably not needed
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


class Neighbourhood:

    def __init__(self, dimensions, max_depth):

        self.dimensions = dimensions
        self.neighbour_table = self.construct_neighbour_table(max_depth)

    def construct_neighbour_table(self, max_depth):
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


class EventHandler:

    def __init__(self):
        pass

    def handle_event(self, event):
        if isinstance(event, BacteriaReplication):
            self.process_bacteria_replication(event)
        elif isinstance(event, RecruitTCell):
            self.process_t_cell_recruitment(event)
        elif isinstance(event, RecruitMacrophage):
            self.process_macrophage_recruitment(event)
        elif isinstance(event, ChemoKillBacteria):
            self.process_chemo_kill_bacteria(event)
        elif isinstance(event, ChemoKillMacrophage):
            self.process_chemo_kill_macrophage(event)

    def process_bacteria_replication(self, event):
        print "BACTERIA REPLICATION"
        self.add_bacteria(event.new_bacteria_address, event.new_metabolism)

    def process_t_cell_recruitment(self, event):
        print "T CELL RECRUITMENT"
        self.add_t_cell(event.t_cell_address)

    def process_macrophage_recruitment(self, event):
        print "MACROPHAGE RECRUITMENT"
        self.add_macrophage(event.macrophage_address, "resting")

    def process_chemo_kill_bacteria(self, event):
        print "CHEMO KILL BACTERIA"
        self.bacteria.remove(event.bacteria_to_kill)
        self.set_attribute_work_grid(event.bacteria_to_kill.address, 'contents', 0.0)

    def process_chemo_kill_macrophage(self, event):
        print "CHEMO KILL MACROPHAGE"
        self.macrophage.remove(event.macrophage_to_kill)
        self.set_attribute_work_grid(event.macrophage_to_kill.address, 'contents', 0.0)


class Automaton(Tile, Neighbourhood, EventHandler):

    def __init__(self, shape, tile_id, attributes, parameters, blood_vessels, fast_bacteria=None, slow_bacteria=None,
                 macrophages=None):
        Tile.__init__(self, shape, attributes)
        Neighbourhood.__init__(self, len(shape), parameters['max_depth'])
        EventHandler.__init__(self)
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
            self.add_bacteria(address, "fast")
        for address in slow_bacteria:
            self.add_bacteria(address, "slow")

    def initialise_macrophages(self, addresses):

        for address in addresses:
            self.add_macrophage(address, "resting")

    def initialise_oxygen_levels(self):

        for address in self.blood_vessels:
            self.set_attribute_grid(address, 'oxygen', self.parameters['initial_oxygen'])

    def update(self):

        self.time += 1

        # ----------------------------
        # CONTINUOUS (Diffusion)
        # ----------------------------

        # Pre-processing (calculating diffusion rates)
        self.diffusion_pre_process()

        if (((self.parameters['chemotherapy_schedule1_start'] / self.parameters['time_step']) <=
                 self.time <
                 (self.parameters['chemotherapy_schedule1_end'] / self.parameters['time_step']))
            or
                (self.parameters['chemotherapy_schedule2_start'] / self.parameters['time_step'] <=
                     self.time)):
            chemo = True

        else:
            chemo = False

        self.max_oxygen_local = 0.0
        self.max_chemotherapy_local = 0.0
        self.max_chemokine_local = 0.0

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
                self.set_attribute_work_grid(address, 'chemotherapy', 0.0)

            # CHEMOKINE
            chemokine_level = self.chemokine(address)
            self.set_attribute_work_grid(address, 'chemokine', chemokine_level)

        # ----------------------------
        # DISCRETE (Agents)
        # ----------------------------

        self.potential_events = []

        # BACTERIA REPLICATION

        for bacteria in self.bacteria:
            bacteria.age += self.parameters['time_step']

            division = False
            if bacteria.metabolism == 'fast':
                max = self.parameters['bacteria_replication_fast_upper']
                min = self.parameters['bacteria_replication_fast_lower']
            else:  # Slow
                max = self.parameters['bacteria_replication_slow_upper']
                min = self.parameters['bacteria_replication_fast_lower']

            replication_time = np.random.randint(min, max)

            if self.time % replication_time == 0:
                division = True

            if division:

                free_neighbours = []

                for depth in range(1,4):
                    if bacteria.neighbourhood == 'mo':
                        neighbours = self.neighbours_moore(bacteria.address, depth, False)
                    else:
                        neighbours = self.neighbours_von_neumann(bacteria.address, depth, False)

                    for neighbour_address in neighbours:
                        neighbour = self.get(neighbour_address)
                        if neighbour is not None and neighbour['contents'] == 0.0 and \
                                neighbour_address not in self.blood_vessels:
                            free_neighbours.append(neighbour_address)

                    if len(free_neighbours) > 0:
                        break

                if len(free_neighbours) == 0:
                    bacteria.resting = True
                    bacteria.age = 0.0
                else: # Free space found
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]

                    if self.address_is_on_grid(neighbour_address):
                        internal = True
                    else:
                        internal = False
                    new_event = BacteriaReplication(neighbour_address, bacteria, internal)
                    self.potential_events.append(new_event)

        # T-CELL RECRUITMENT
        if self.number_of_bacteria_global >= self.parameters['bacteria_threshold_for_t_cells']:
            # Each blood vessel
            for bv_address in self.blood_vessels:
                r = np.random.randint(1, 100)
                if r <= self.parameters['t_cell_recruitment_probability']:
                    neighbours = self.neighbours_von_neumann(bv_address, 1)
                    # Get neighbours which are empty and have a high enough chemokine level
                    free_neighbours = []
                    for n in neighbours:
                        if self.get(n) is not None and \
                                self.get_attribute(n, 'blood_vessel') == 0.0 and \
                                self.get_attribute(n, 'contents') == 0.0 and \
                                self.chemokine_scale(n) > self.parameters['chemokine_scale_for_t_cell_recruitment']:
                            free_neighbours.append(n)
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

        # MACROPHAGE RECRUITMENT
        for bv_address in self.blood_vessels:
            r = np.random.randint(1, 100)
            if r <= self.parameters['macrophage_recruitment_probability']:
                neighbours = self.neighbours_von_neumann(bv_address, 1)
                free_neighbours = []
                for n in neighbours:
                    if self.get(n) is not None and \
                                self.get_attribute(n, 'blood_vessel') == 0.0 and \
                                self.get_attribute(n, 'contents') == 0.0 and \
                                self.chemokine_scale(n) > self.parameters['chemokine_scale_for_macrophage_recruitment']:
                        free_neighbours.append(n)

                if len(free_neighbours) > 0:
                    # Pick one of the neighbours
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                    if self.address_is_on_grid(neighbour_address):
                        internal = True
                    else:
                        internal = False
                    new_event = RecruitMacrophage(neighbour_address, internal)
                    self.potential_events.append(new_event)

        # CHEMOTHERAPY KILLING BACTERIA
        for b in self.bacteria:
            chemo_scale = self.chemotherapy_scale(b.address)
            if (b.metabolism == 'fast' and chemo_scale > self.parameters['chemotherapy_scale_for_kill_fast_bacteria']) \
                or \
                (b.metabolism == 'slow' and chemo_scale > self.parameters['chemotherapy_scale_for_kill_slow_bacteria']):
                new_event = ChemoKillBacteria(b)
                self.potential_events.append(new_event)

        # CHEMOTHERAPY KILLING MACROPHAGES
        for m in self.macrophages:
            chemo_scale = self.chemotherapy_scale(m.address)
            if ((m.state == 'infected' or m.state == 'chronically_infected') and chemo_scale >
                    self.parameters['chemotherapy_scale_for_kill_macrophage']):
                new_event = ChemoKillMacrophage(m)
                self.potential_events.append(new_event)

        # TODO - T-cell movement & killing

        # TODO - Macrophage movement & death

        # Reorder events
        self.reorder_events()

    def reorder_events(self):
        # TODO - other methods - currently just random
        np.random.shuffle(self.potential_events)

    def process_events(self, events):

        for event in events:
            self.handle_event(event)

        self.persist_agents()

        self.swap_grids()

    def diffusion_pre_process(self):

        for location in range(self.size):
            address = self.location_to_address(location)

            # On grid
            oxygen_diffusion = self.parameters['oxygen_diffusion']
            chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

            # Check if there is specified amount of caseum within specified distance of cell
            neighbours = self.neighbours_moore(address, int(self.parameters['caseum_distance']))
            caseum_count = 0
            for neighbour_address in neighbours:
                cell = self.get(neighbour_address)
                if cell is not None and cell['contents'] == 'caseum':
                    caseum_count += 1
                    # Once the caseum threshold is reached
                    if caseum_count == self.parameters['caseum_threshold']:
                        # Decrease the diffusion level at the cell
                        oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                        chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                        # Exit the loop
                        break

            # Need to set the values on the current grid
            self.set_attribute_grid(address, 'oxygen_diffusion_rate', oxygen_diffusion)
            self.set_attribute_grid(address, 'chemotherapy_diffusion_rate', chemotherapy_diffusion)

        # Set diffusion rates on halo
        for halo_address in self.halo_depth1:
            if self.get(halo_address) is not None:

                # On grid
                oxygen_diffusion = self.parameters['oxygen_diffusion']
                chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']

                neighbours = self.neighbours_moore(halo_address, int(self.parameters['caseum_distance']))
                caseum_count = 0
                for neighbour_address in neighbours:
                    cell = self.get(neighbour_address)
                    if cell is not None and cell['contents'] == 'caseum':
                        caseum_count += 1
                        # Once the caseum threshold is reached
                        if caseum_count == self.parameters['caseum_threshold']:
                            # Decrease the diffusion level at the cell
                            oxygen_diffusion /= self.parameters['oxygen_diffusion_caseum_reduction']
                            chemotherapy_diffusion /= self.parameters['chemotherapy_diffusion_caseum_reduction']
                            # Exit the loop
                            break

                # Need to set the values on the current grid
                index = self.halo_addresses.index(halo_address)
                self.halo_cells[index]['oxygen_diffusion_rate'] = oxygen_diffusion
                self.halo_cells[index]['chemotherapy_diffusion_rate'] = chemotherapy_diffusion

    def oxygen(self, address):
        # Get the current cell values
        cell = self.get(address)

        # Get diffusion value for cell
        cell_diffusion = cell['oxygen_diffusion_rate']

        # Initialise expression
        expression = 0

        # Get immediate von neumann neighbours
        neighbour_addresses = self.neighbours_von_neumann(address,1)
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
        if isinstance(cell['contents'], Bacteria):
            expression -= self.parameters['oxygen_uptake_from_bacteria'] * cell['oxygen']

        # Calculate new level
        new_oxygen = cell['oxygen'] + self.parameters['time_step'] * expression

        # Overwrite the maximum oxygen value if larger
        self.max_oxygen_local = max(self.max_oxygen_local, new_oxygen)

        return new_oxygen

    def chemotherapy(self, address):
        # Get the current cell values
        cell = self.get(address)

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
        cell = self.get(address)

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
        if isinstance(cell['contents'], Bacteria):
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

    def add_bacteria(self, address, metabolism):
        new_bacteria = Bacteria(address, metabolism)
        self.bacteria.append(new_bacteria)
        self.set_attribute_work_grid(address, 'contents', new_bacteria)

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


# ------------------------------
# AGENTS
# ------------------------------
class Agent:

    def __init__(self, address):
        self.address = address
        self.age = 0


class Bacteria(Agent):

    def __init__(self, address, metabolism):
        self.metabolism = metabolism
        self.neighbourhood = 'mo'
        self.resting = False
        Agent.__init__(self, address)


class Macrophage(Agent):

    def __init__(self, address, state):
        self.state = state
        Agent.__init__(self, address)


class TCell(Agent):

    def __init__(self, address):
        Agent.__init__(self, address)
        pass


# ------------------------------
# EVENTS
# ------------------------------
class Event:

    def __init__(self, addresses_affected, internal):
        self.addresses_affected = addresses_affected
        self.internal = internal


class BacteriaReplication(Event):

    def __init__(self, address, bacteria, internal):
        Event.__init__(self, [address], internal)
        self.new_bacteria_address = address
        self.new_metabolism = bacteria.metabolism
        self.original_bacteria = bacteria

    def clone(self, new_addresses):
        return BacteriaReplication(new_addresses[0], self.original_bacteria, self.new_metabolism)


class RecruitTCell(Event):

    def __init__(self, address, internal):
        self.t_cell_address = address
        Event.__init__(self,[address], internal)

    def clone(self, new_addresses):
        return RecruitTCell(new_addresses[0], self.internal)


class RecruitMacrophage(Event):

    def __init__(self, address, internal):
        self.macrophage_address = address
        Event.__init__(self, [address], internal)

    def clone(self, new_addresses):
        return RecruitMacrophage(new_addresses[0], self.internal)


class ChemoKillBacteria(Event):

    def __init__(self, bacteria_to_kill):
        self.bacteria_to_kill = bacteria_to_kill
        # Chemo killing is always internal
        Event.__init__(self, [bacteria_to_kill.address], True)


class ChemoKillMacrophage(Event):

    def __init__(self, macrophage_to_kill):
        self.macrophage_to_kill = macrophage_to_kill
        # Chemo killing is always internal
        Event.__init__(self, [macrophage_to_kill.address], True)