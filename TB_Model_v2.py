import numpy as np
from collections import Counter
import itertools
import math


class Tile:

    def __init__(self, shape, attributes, dtype, max_depth):
        pass
        self.attributes = attributes
        self.shape = shape
        self.max_depth = max_depth

        dtype.append(('neighbours_mo', np.object))
        dtype.append(('neighbours_vn', np.object))

        self.grid = np.zeros(shape, dtype=dtype)
        self.work_grid = None

        self.calculate_neighbours_2D()

    def create_work_grid(self):
        self.work_grid = self.grid.copy()

    def swap_grids(self):
        self.grid, self.work_grid = self.work_grid, self.grid

    def set_attribute_grid(self,address,attribute,value):
        self.grid[address][attribute] = value

    def set_attribute_work_grid(self,address,attribute,value):
        self.work_grid[address][attribute] = value

    def calculate_neighbours_2D(self):

        # Initialise empty dictionaries
        self.moore_relative = dict()
        self.von_neumann_relative = dict()
        # Add an entry for each depth
        for d in range(1, self.max_depth + 1):
            self.von_neumann_relative[d] = []

        for depth in range(1, self.max_depth + 1):
            # Get truth table values (e.g. depth 2 gives [-2,-1,0,1,2] for range_)
            range_ = range(-depth, depth + 1)
            # Use product to find all combinations for given depth and number of dimensions
            row = list(itertools.product(range_, repeat=2))
            # Remove the 0 entry e.g. (0,0) for 2 dimensions
            row.remove((0,) * 2)

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

        it = np.nditer(self.grid, flags=['multi_index','refs_ok'], op_flags=['writeonly'])
        while not it.finished:

            x = it.multi_index[0]
            y = it.multi_index[1]

            vn = dict()
            mo = dict()

            for d in range(1, self.max_depth+1):

                vn[d] = []
                # VN
                vn_row = self.von_neumann_relative[d]
                for address in vn_row:
                    relative_neighbour = (x + address[0], y + address[1])
                    if self.address_is_on_grid(relative_neighbour):
                        vn[d].append(relative_neighbour)

                # MO
                mo[d] = []
                # VN
                mo_row = self.moore_relative[d]
                for address in mo_row:
                    relative_neighbour = (x + address[0], y + address[1])
                    if self.address_is_on_grid(relative_neighbour):
                        mo[d].append(relative_neighbour)

            it[0]['neighbours_vn'] = vn
            it[0]['neighbours_mo'] = mo

            it.iternext()

    def address_is_on_grid(self, address):
        for a in range(len(address)):
            if address[a] < 0 or address[a] > self.shape[a]:
                return False
        return True


class Automaton(Tile):

    def __init__(self, shape, parameters, blood_vessels, fast_bacteria=None, slow_bacteria=None,
                 macrophages=None, output_location=''):

        attributes = ['oxygen', 'chemotherapy', 'chemokine', 'contents', 'blood_vessel', 'oxygen_diffusion_rate',
                      'chemotherapy_diffusion_rate']

        dtype = []
        dtype.append(('oxygen', np.double))
        dtype.append(('chemotherapy', np.double))
        dtype.append(('chemokine', np.double))
        dtype.append(('blood_vessel', np.double))
        dtype.append(('oxygen_diffusion_rate', np.double))
        dtype.append(('chemotherapy_diffusion_rate', np.double))
        dtype.append(('contents', np.object))

        self.parameters = parameters
        Tile.__init__(self, shape, attributes, dtype, self.parameters['max_depth'])

        self.time = self.parameters['initial_time']

        # LISTS
        self.bacteria = []
        self.macrophages = []
        self.t_cells = []
        self.caseum_addresses = []
        self.blood_vessel_addresses = []
        self.potential_events = []

        self.max_oxygen = 0.0
        self.max_chemotherapy = 0.0
        self.max_chemokine = 0.0

        self.initialise(blood_vessels, fast_bacteria, slow_bacteria, macrophages)
        self.create_work_grid()

        self.chemo_schedule1_start = np.random.randint(self.parameters['chemotherapy_schedule1_start_lower'],
                                                       self.parameters['chemotherapy_schedule1_start_upper'])

        # Set up output file paths
        if output_location != '':
            output_location += '/'
        self.totalcell_test_file_path = output_location + 'totalcell_test.txt'
        self.contents_file_path = output_location + 'data_test.txt'
        self.oxygen_file_path = output_location + 'oxygen_test.txt'
        self.chemotherapy_file_path = output_location + 'chemo1.txt'
        self.chemokine_file_path = output_location + 'ckine.txt'
        self.type1_file_path = output_location + 'Type1.txt'
        self.type1_r_file_path = output_location + 'Type1_R.txt'
        self.type2_file_path = output_location + 'Type2.txt'
        self.type2_r_file_path = output_location + 'Type2_R.txt'
        self.type3_file_path = output_location + 'Type3.txt'
        self.activemac_file_path = output_location + 'activemac.txt'
        self.restingmac_file_path = output_location + 'restingmac.txt'
        self.infectedmac_file_path = output_location + 'infectedmac.txt'
        self.chroninfectedmac_file_path = output_location + 'chroninfectedmac.txt'
        self.caseation_file_path = output_location + 'caseation.txt'
        self.total_file_path = output_location + 'Total.txt'
        self.intra_bac_file_path = output_location + 'intra_bac.txt'

        self.type1_file = open(self.type1_file_path, 'w')
        self.type1_r_file = open(self.type1_r_file_path, 'w')
        self.type2_file = open(self.type2_file_path, 'w')
        self.type2_r_file = open(self.type2_r_file_path, 'w')
        self.type3_file = open(self.type3_file_path, 'w')
        self.activemac_file = open(self.activemac_file_path, 'w')
        self.restingmac_file = open(self.restingmac_file_path, 'w')
        self.infectedmac_file = open(self.infectedmac_file_path, 'w')
        self.chroninfectedmac_file = open(self.chroninfectedmac_file_path, 'w')
        self.caseation_file = open(self.caseation_file_path, 'w')
        self.total_file = open(self.total_file_path, 'w')
        self.total_cell_test_file = open(self.totalcell_test_file_path, 'w')
        self.intra_bac_file = open(self.intra_bac_file_path, 'w')

        self.contents_file = open(self.contents_file_path, 'w')
        self.oxygen_file = open(self.oxygen_file_path, 'w')
        self.chemotherapy_file = open(self.chemotherapy_file_path, 'w')
        self.chemokine_file = open(self.chemokine_file_path, 'w')

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

    def initialise(self, blood_vessels, fast_bacteria, slow_bacteria, macrophages):

        for bv_address in blood_vessels:
            self.set_attribute_grid(bv_address, 'blood_vessel', self.parameters['blood_vessel_value'])
            self.blood_vessel_addresses.append(bv_address)
            self.set_attribute_grid(bv_address, 'oxygen', self.parameters['initial_oxygen']
                                    * self.parameters['blood_vessel_value'])
        self.max_oxygen = self.parameters['initial_oxygen'] * self.parameters['blood_vessel_value']
        for fb_address in fast_bacteria:
            fast_bacterium = Bacterium(fb_address, 'fast')
            self.bacteria.append(fast_bacterium)
            self.set_attribute_grid(fb_address, 'contents', fast_bacterium)
        for sb_address in slow_bacteria:
            slow_bacterium = Bacterium(sb_address, 'slow')
            self.bacteria.append(slow_bacterium)
            self.set_attribute_grid(sb_address, 'contents', slow_bacterium)
        for m_address in macrophages:
            macrophage = Macrophage(m_address, 'resting')
            self.macrophages.append(macrophage)
            self.set_attribute_grid(m_address, 'contents', macrophage)

        ita = np.nditer(self.grid, flags=['multi_index', 'refs_ok'], op_flags=['writeonly'])
        while not ita.finished:
            ita[0]['oxygen_diffusion_rate'] = self.parameters['oxygen_diffusion']
            ita[0]['chemotherapy_diffusion_rate'] = self.parameters['chemotherapy_diffusion']
            ita.iternext()

    def run(self):

        while self.time * self.parameters['time_step'] <= self.parameters['time_limit']:

            self.time += 1.0

            self.update()

            self.resolve_conflict()

            self.process_events()

            
    def update(self):

        # CONTINUOUS
        self.diffusion_pre_process()

        chemo = (self.chemo_schedule1_start / self.parameters['time_step']) <= self.time < \
                (self.parameters['chemotherapy_schedule1_end'] / self.parameters['time_step']) or \
                self.parameters['chemotherapy_schedule2_start'] / self.parameters['time_step'] <= self.time

        self.diffusion(chemo)


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
        caseum_addresses = list(self.caseum_addresses)

        # Process all caseum addresses
        for address in caseum_addresses:
            # Check the nighbourhood of the cell, up to the pre-set depth
            for depth in range(1, int(self.parameters['caseum_distance_to_reduce_diffusion']+1)):
                # Record each affected neighbour in the list (can be duplicates in list)
                neighbours = self.grid[address]['neighbours_moore'][depth]
                for neighbour in neighbours:
                    affected_addresses.append(neighbour)

        # Count how many times each address appears in the affected list
        counted = Counter(affected_addresses)

        # Loop through every address
        for address in counted:
            grid_cell = self.grid[address]
            # Get initial diffusion rates
            oxygen_diffusion = self.parameters['oxygen_diffusion']
            chemotherapy_diffusion = self.parameters['chemotherapy_diffusion']
            # If the address has been affected and the number of affectations exceeds the threshold
            if counted[address] >= self.parameters['caseum_threshold_to_reduce_diffusion']:
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

    def diffusion(self, chemo):
        above = np.roll(self.grid, shift=1, axis=0)
        below = np.roll(self.grid, shift=-1, axis=0)
        left = np.roll(self.grid, shift=1, axis=1)
        right = np.roll(self.grid, shift=-1, axis=1)

        # MAIN GRID - excluding edges
        # OXYGEN
        new_oxygen = self.grid['oxygen'] + self.parameters['time_step'] \
                         * (((((self.grid['oxygen_diffusion_rate']+ below['oxygen_diffusion_rate']) / 2)
                         * (below['oxygen'] - self.grid['oxygen']))
                         - (((self.grid['oxygen_diffusion_rate'] + above['oxygen_diffusion_rate']) / 2)
                         * (self.grid['oxygen'] - above['oxygen'])))
                         / self.parameters['spatial_step'] ** 2
                         + ((((self.grid['oxygen_diffusion_rate'] + right['oxygen_diffusion_rate']) / 2)
                         * (right['oxygen'] - self.grid['oxygen']))
                         - (((self.grid['oxygen_diffusion_rate']+ left['oxygen_diffusion_rate']) / 2)
                         * (self.grid['oxygen'] - left['oxygen'])))
                         / self.parameters['spatial_step'] ** 2
                         + self.parameters['oxygen_from_source'] * self.grid['blood_vessel']
                         - self.parameters['oxygen_uptake_from_bacteria'] * self.grid['oxygen']
                         * isinstance(self.grid['contents'],Bacterium))

        self.work_grid['oxygen'][1:self.work_grid.shape[0]-1] = new_oxygen[1:self.work_grid.shape[0]-1]
        
        # CHEMOTHERAPY
        if chemo:
            new_chemotherapy = self.grid['chemotherapy'] + self.parameters['time_step'] \
                                   * (((((self.grid['chemotherapy_diffusion_rate'] + below['chemotherapy_diffusion_rate']) / 2)
                                   * (below['chemotherapy'] - self.grid['chemotherapy']))
                                   - (((self.grid['chemotherapy_diffusion_rate'] + above['chemotherapy_diffusion_rate']) / 2)
                                   * (self.grid['chemotherapy'] - above['chemotherapy'])))
                                   / self.parameters['spatial_step'] ** 2
                                   + ((((self.grid['chemotherapy_diffusion_rate'] + right['chemotherapy_diffusion_rate']) / 2)
                                   * (right['chemotherapy'] - self.grid['chemotherapy']))
                                   - (((self.grid['chemotherapy_diffusion_rate'] + left['chemotherapy_diffusion_rate']) / 2)
                                   * (self.grid['chemotherapy'] - left['chemotherapy'])))
                                   / self.parameters['spatial_step'] ** 2
                                   + self.parameters['chemotherapy_from_source'] * self.grid['blood_vessel']
                                   - self.parameters['chemotherapy_decay'] * self.grid['chemotherapy'])
            self.work_grid['chemotherapy'][1:self.work_grid.shape[0] - 1] = new_chemotherapy[1:self.work_grid.shape[0] - 1]
        else:
            self.work_grid['chemotherapy'] = np.zeros(self.grid.shape)

        # CHEMOKINE
        new_chemokine = self.grid['chemokine'] + self.parameters['time_step'] \
                            * ((((self.parameters['chemokine_diffusion'])
                            * (below['chemokine'] - self.grid['chemokine']))
                            - ((self.parameters['chemokine_diffusion'])
                            * (self.grid['chemokine'] - above['chemokine'])))
                            / self.parameters['spatial_step'] ** 2
                            + (((self.parameters['chemokine_diffusion'])
                            * (right['chemokine'] - self.grid['chemokine']))
                            - ((self.parameters['chemokine_diffusion'])
                            * (self.grid['chemokine'] - left['chemokine'])))
                            / self.parameters['spatial_step'] ** 2
                            + (self.parameters['chemokine_from_bacteria'] * isinstance(self.grid['contents'], Bacterium))
                            + (self.parameters['chemokine_from_macrophage']
                            * (isinstance(self.grid['contents'], Macrophage) and self.grid['contents'].state !='resting'))
                            - self.parameters['chemokine_decay'] * self.grid['chemokine'])

        self.work_grid['chemokine'][1:self.work_grid.shape[0] - 1] = new_chemokine[1:self.work_grid.shape[0] - 1]

        # Function to diffuse into the corner cells
        def corner_diffusion(address, neighbour1, neighbour2):
            new_oxygen = self.grid[address]['oxygen'] + self.parameters['time_step'] \
                         * (self.grid[address]['oxygen_diffusion_rate']
                         * (self.grid[neighbour1]['oxygen'] - 2 * self.grid[address]['oxygen'] + self.grid[neighbour1]['oxygen'])
                         / (self.parameters['spatial_step'] ** 2)
                         + self.grid[address]['oxygen_diffusion_rate']
                         * (self.grid[neighbour2]['oxygen'] - 2 * self.grid[address]['oxygen'] + self.grid[neighbour2]['oxygen'])
                         / (self.parameters['spatial_step'] ** 2)
                         + self.parameters['oxygen_from_source'] * self.grid[address]['blood_vessel']
                         - (self.parameters['oxygen_uptake_from_bacteria'] *self.grid[address]['oxygen']
                         * isinstance(self.grid[address]['contents'], Bacterium)))
            self.work_grid[address]['oxygen'] = new_oxygen
            if chemo:
                new_chemotherapy = self.grid[address]['chemotherapy'] + self.parameters['time_step'] \
                                    * (self.grid[address]['chemotherapy_diffusion_rate']
                                    * (self.grid[neighbour1]['chemotherapy'] - 2 * self.grid[address]['chemotherapy']
                                    + self.grid[neighbour1]['chemotherapy'])
                                    / (self.parameters['spatial_step'] ** 2)
                                    + self.grid[address]['chemotherapy_diffusion_rate']
                                    * (self.grid[neighbour2]['chemotherapy'] - 2 * self.grid[address]['chemotherapy']
                                    + self.grid[neighbour2]['chemotherapy'])
                                    / (self.parameters['spatial_step'] ** 2)
                                    + self.parameters['chemotherapy_from_source'] * self.grid[address]['blood_vessel']
                                    - self.parameters['chemotherapy_decay'] * self.grid[address]['chemotherapy'])
                self.work_grid[address]['chemotherapy'] = new_chemotherapy
            new_chemokine = self.grid[address]['oxygen'] + self.parameters['time_step'] \
                            * (self.grid[address]['oxygen_diffusion_rate']
                            * (self.grid[neighbour1]['oxygen'] - 2 * self.grid[address]['oxygen']
                            + self.grid[neighbour1]['oxygen'])
                            / (self.parameters['spatial_step'] ** 2)
                            + self.grid[address]['oxygen_diffusion_rate']
                            * (self.grid[neighbour2]['oxygen'] - 2 * self.grid[address]['oxygen'] + self.grid[neighbour2]['oxygen'])
                            / (self.parameters['spatial_step'] ** 2)
                            + (self.parameters['chemokine_from_bacteria'] * isinstance(self.grid[address]['contents'], Bacterium))
                            + (self.parameters['chemokine_from_macrophage'] * (isinstance(self.grid['contents'], Macrophage) and self.grid[address]['contents'].state != 'resting'))
                            - self.parameters['chemokine_decay'] * self.grid[address]['chemokine'])
            self.work_grid[address]['chemokine'] = new_chemokine

        corner_diffusion((0, 0), (1, 0), (0, 1))
        corner_diffusion((0, self.shape[1]-1), (1, self.shape[1]-1), (0, self.shape[1]-2))
        corner_diffusion((self.shape[0]-1, 0), (self.shape[0]-2, 0), (self.shape[0]-1, 1))
        corner_diffusion((self.shape[0]-1, self.shape[1]-1), (self.shape[0]-2, self.shape[1]-1), (self.shape[0]-1, self.shape[1]-2))

        # TOP ROW
        def edge_diffusion(row, values1, values2, values3):
            """
            
            :param row: 
            :param values1: Direction which does not have corresponding alternative (i.e. Below for top row) 
            :param values2: 
            :param values3: 
            :return: 
            """
            new_oxygen = row['oxygen'] + self.parameters['time_step'] \
                    * (row['oxygen_diffusion_rate'] * (values1['oxygen'] - 2 * row['oxygen'] + values1['oxygen']) 
                    / (self.parameters['spatial_step']**2) 
                    + row['oxygen_diffusion_rate'] * (values2['oxygen'] - 2 * row['oxygen'] + values3['oxygen']) 
                    / (self.parameters['spatial_step']**2)
                    + self.parameters['oxygen_from_source'] * row['blood_vessel']
                    - self.parameters['oxygen_uptake_from_bacteria'] * row['oxygen']
                    * isinstance(row['contents'], Bacterium))
            if chemo:
                new_chemotherapy = row['chemotherapy'] + self.parameters['time_step'] \
                                    * (row['chemotherapy_diffusion_rate'] * (
                                    values1['chemotherapy'] - 2 * row['chemotherapy'] + values1['chemotherapy'])
                                    / (self.parameters['spatial_step'] ** 2)
                                    + row['chemotherapy_diffusion_rate'] 
                                    * (values2['chemotherapy'] - 2 * row['chemotherapy'] + values3['chemotherapy'])
                                    / (self.parameters['spatial_step'] ** 2)
                                    + self.parameters['chemotherapy_from_source'] * row['blood_vessel']
                                    - self.parameters['chemotherapy_decay'] * row['chemotherapy'])
            else:
                new_chemotherapy = np.zeros(row.shape)

            new_chemokine = row['chemokine'] + self.parameters['time_step'] \
                            * (self.parameters['chemokine_diffusion']
                            * (values1['chemokine'] - 2 * row['chemokine'] + values1['chemokine'])
                            / (self.parameters['spatial_step'] ** 2)
                            + self.parameters['chemokine_diffusion']
                            * (values2['chemokine'] - 2 * row['chemokine'] + values3['chemokine'])
                            / (self.parameters['spatial_step'] ** 2)
                            + (self.parameters['chemokine_from_bacteria'] * isinstance(row['contents'], Bacterium))
                            + (self.parameters['chemokine_from_macrophage'] * (isinstance(row['contents'], Macrophage) and row['contents'].state != 'resting'))
                            - self.parameters['chemokine_decay'] * row['chemokine'])

            return new_oxygen, new_chemotherapy, new_chemokine

        self.work_grid[0, 1:self.shape[1] - 1]['oxygen'], self.work_grid[0, 1:self.shape[1] - 1]['chemotherapy'], \
        self.work_grid[0, 1:self.shape[1] - 1]['chemokine'] \
            = edge_diffusion(self.grid[0, 1:self.shape[1] - 1],
                             below[0, 1:self.shape[1] - 1],
                             left[0, 1:self.shape[1] - 1],
                             right[0, 1:self.shape[1] - 1])

        # BOTTOM ROW
        self.work_grid[self.shape[0]-1, 1:self.shape[1] - 1]['oxygen'], self.work_grid[self.shape[0]-1, 1:self.shape[1] - 1]['chemotherapy'], \
        self.work_grid[self.shape[0]-1, 1:self.shape[1] - 1]['chemokine'] \
            = edge_diffusion(self.grid[self.shape[0]-1, 1:self.shape[1] - 1],
                             above[self.shape[0]-1, 1:self.shape[1] - 1],
                             left[self.shape[0]-1, 1:self.shape[1] - 1],
                             right[self.shape[0]-1, 1:self.shape[1] - 1])

        # LEFT COLUMN
        self.work_grid[1:self.shape[0]-1, 0]['oxygen'], self.work_grid[1:self.shape[0]-1, 0]['chemotherapy'], \
        self.work_grid[1:self.shape[0]-1, 0]['chemokine'] \
            = edge_diffusion(self.grid[1:self.shape[0]-1, 0],
                             right[1:self.shape[0]-1, 0],
                             above[1:self.shape[0]-1, 0],
                             below[1:self.shape[0]-1, 0])

        # RIGHT COLUMN
        self.work_grid[1:self.shape[0] - 1, self.shape[1] - 1]['oxygen'], self.work_grid[1:self.shape[0] - 1, self.shape[1] - 1]['chemotherapy'], \
        self.work_grid[1:self.shape[0] - 1, self.shape[1] - 1]['chemokine'] \
            = edge_diffusion(self.grid[1:self.shape[0] - 1, self.shape[1] - 1],
                             left[1:self.shape[0] - 1, self.shape[1] - 1],
                             above[1:self.shape[0] - 1, self.shape[1] - 1],
                             below[1:self.shape[0] - 1, self.shape[1] - 1])


class Agent:

    def __init__(self, address):
        self.address = address


class Bacterium(Agent):

    def __init__(self, address, metabolism):
        Agent.__init__(self, address)
        self.metabolism = metabolism


class Macrophage(Agent):

    def __init__(self, address, state):
        Agent.__init__(self, address)
        self.state = state


