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

    def get_total_bacteria(self):
        return len(self.bacteria) + sum([m.intracellular_bacteria for m in self.macrophages])

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
            self.max_oxygen = self.grid['oxygen'].max()
            self.max_chemotherapy = self.grid['chemotherapy'].max()
            self.max_chemokine = self.grid['chemokine'].max()

            self.update()

            self.resolve_conflicts()

            self.process_events()

    def update(self):

        # ----------------------------
        # CONTINUOUS (Diffusion)
        # ----------------------------
        self.diffusion_pre_process()

        chemo = (self.chemo_schedule1_start / self.parameters['time_step']) <= self.time < \
                (self.parameters['chemotherapy_schedule1_end'] / self.parameters['time_step']) or \
                self.parameters['chemotherapy_schedule2_start'] / self.parameters['time_step'] <= self.time

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

            # If the time is sufficient enough, bacteria can replicate
            if self.time % replication_time == 0:

                # Look for free neighbours
                free_neighbours = []
                for depth in range(1, self.parameters['max_depth']+1):
                    # Pull the neighbours from the appropriate neighbourhood
                    if bacterium.division_neighbourhood == 'mo':
                        neighbours = self.grid[bacterium.address]['neighbour_mo'][depth]
                    else:
                        neighbours = self.grid[bacterium.address]['neighbour_vn'][depth]
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
                    event_details = dict(type='resting',value=True)
                    new_event = Event('BacteriumStateChange', [bacterium.address], [bacterium.address], event_details)
                    self.potential_events.append(new_event)
                else:  # Free space found
                    # Pick a free neighbour at random
                    neighbour_address = free_neighbours[np.random.randint(len(free_neighbours))]
                    event_details = dict(from_address=bacterium.address, to_address=neighbour_address,
                                         metabolism=bacterium.metabolsim)
                    new_event = Event('BacteriumReplication', [bacterium.address, neighbour_address],
                                      [bacterium.address, neighbour_address], event_details)
                    self.potential_events.append(new_event)

    def t_cell_recruitment(self):
        """
        Once bacteria over entire system reach a threshold, t-cells enter the system. Creates an event to add a t-cell
        to a cell next to a blood vessel
        :return:
        """
        # When global amount of bacteria exceeds threshold
        if self.get_total_bacteria() >= self.parameters['bacteria_threshold_for_t_cells']:
            # Each blood vessel
            for bv_address in self.blood_vessel_addresses:
                # Generate event if probability according to parameters
                r = np.random.randint(1, 101)
                if r <= self.parameters['t_cell_recruitment_probability']:
                    # Get von Neumann neighbours
                    neighbours = self.grid[bv_address]['neighbours_vn'][1]
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
                        # Create event (no extra details needed outside the address)
                        new_event = Event('RecruitTCell', [neighbour_address], [neighbour_address], None)
                        self.potential_events.append(new_event)

    def macrophage_recruitment(self):
        """
        Each step for each source vessel, there is a probability that macrophage will be recruited
        :return:
        """
        if self.get_total_bacteria() >= self.parameters['bacteria_threshold_for_macrophage_recruitment']:
            chemokine_threshold = self.parameters['chemokine_scale_for_macrophage_recruitment_above_threshold']
        else:
            chemokine_threshold = self.parameters['chemokine_scale_for_macrophage_recruitment_below_threshold']

        # Loop through each blood vessel
        for bv_address in self.blood_vessel_addresses:
            # Generate event with probability based on parameters
            r = np.random.randint(1, 101)
            if r <= self.parameters['macrophage_recruitment_probability']:
                # Get neighbours, then reduce to those that are free and have sufficient chemokine scale
                neighbours = self.grid[bv_address]['neighbours_vn'][1]
                free_neighbours = []
                for neighbour_address in neighbours:
                    neighbour = self.grid[neighbour_address]
                    if neighbour is not None and neighbour['blood_vessel'] == 0.0 and neighbour['contents'] == 0.0 and \
                            self.chemokine_scale(neighbour_address) > chemokine_threshold:
                        free_neighbours.append(neighbour_address)

                if len(free_neighbours) > 0:
                    # Pick one of the neighbours
                    chosen_neighbour = free_neighbours[np.random.randint(len(free_neighbours))]
                    # Create event (no details needed outside of the address)
                    new_event = Event('RecruitMacrophage', [chosen_neighbour], [chosen_neighbour], None)
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
                    or (bacterium.metabolism == 'slow' and chemo_scale >
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
        pass

    def oxygen_scale(self, address):
        """
        Oxygen level at cell as % of maximum global oxygen level
        :param address:
        :return:
        """
        if self.max_oxygen == 0.0:
            return 0.0
        else:
            return (self.grid[address]['oxygen'] / self.max_oxygen) * 100

    def chemotherapy_scale(self, address):
        """
        Maximum chemotherapy level at cell as % of global maximum chemotherapy
        :param address:
        :return:
        """
        if self.max_chemotherapy == 0.0:
            return 0.0
        else:
            return (self.grid[address]['chemotherapy'] / self.max_chemotherapy) * 100

    def chemokine_scale(self, address):
        """
        Chemokine level at cell as % of global maximum chemokine level
        :param address:
        :return:
        """
        if self.max_chemokine == 0.0:
            return 0.0
        else:
            return (self.grid[address]['chemokine'] / self.max_chemokine) * 100.0


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
        self.intracellular_bacteria = 0


class Event:

    def __init__(self, event_type, dependant_addresses, impacted_addresses, event_details):
        self.event_type = event_type
        self.dependant_addresses = dependant_addresses
        self.impacted_addresses = impacted_addresses
        self.allowed_impacted_addresses = []
        self.event_details = event_details


