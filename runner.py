import TB_Model
import ConfigParser


def run_single(topology, time_limit):
    pass


def run_many_serial(topology, time_limit):

    number_tiles = len(topology.automata)

    # HALOS
    dz_addresses = []
    for tile_id in range(number_tiles):
        dz_addresses.append(topology.automata[tile_id].danger_zone_addresses)

    values = []
    halo_addresses = topology.external_addresses_required

    for t in range(time_limit):

        print "TIME-STEP:", t

        for tile_id in range(number_tiles):
            # ---------------- BACK -----------------------------
            values.append(topology.automata[tile_id].get_danger_zone())
            # ---------------- BACK -----------------------------

        halos = construct_halos(topology, dz_addresses, values, halo_addresses)
        max_oxygen = 0.0
        max_chemotherapy = 0.0
        max_chemokine = 0.0
        for automaton in topology.automata:
            max_oxygen = max(max_oxygen, automaton.max_oxygen_local)
            max_chemotherapy = max(max_chemotherapy, automaton.max_chemotherapy_local)
            max_chemokine = max(max_chemokine, automaton.max_chemokine_local)

        for automaton in topology.automata:
            print "Running automata: ", automaton.tile_id

            # UPDATE
            # ------------------ OUT ----------------------------
            automaton.set_halo(halos[automaton.tile_id])
            automaton.set_max_oxygen_global(max_oxygen)
            automaton.set_max_chemotherapy_global(max_chemotherapy)
            automaton.set_max_chemokine_global(max_chemokine)
            # ------------------ OUT ----------------------------

            automaton.update()

            print automaton.potential_events


def run_many_parallel(topology, time_limit):
    pass


def construct_halos(topology, danger_zone_addresses, danger_zone_values, halo_addresses):

    number_tiles = len(danger_zone_addresses)

    halos = []

    for tile_id in range(number_tiles):
        halos.append([])

    for tile_id in range(number_tiles):
        for local_address in halo_addresses:
            global_address = topology.local_to_global(tile_id, local_address)

            new_tile_id, new_local_address = topology.global_to_local(global_address)
            if new_tile_id is not None:
                index = danger_zone_addresses[new_tile_id].index(new_local_address)
                halos[tile_id].append(danger_zone_values[new_tile_id][index])
            else:
                halos[tile_id].append(None)

    return halos


def initialise(config):

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")

    blood_vessel_addresses = []
    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", "blood_vessels_hard_code").split("/")
        for b in bv_list:
            blood_vessel_addresses.append([int(c) for c in b.split(",")])
    # TODO - random initialise

    # FAST BACTERIA
    bacteria_method = config.get("InitialiseSection","bacteria")

    if bacteria_method == 'hard_code':
        fast_list = config.get("InitialiseSection", "bacteria_fast_hard_code").split("/")
        fast_addresses = []
        for a in fast_list:
            address = [int(c) for c in a.split(",")]
            if address not in blood_vessel_addresses:
                fast_addresses.append(address)
            else:
                # TODO - avoid conflict
                pass
        slow_list = config.get("InitialiseSection", "bacteria_slow_hard_code").split("/")
        slow_addresses = []
        for a in slow_list:
            address = [int(c) for c in a.split(",")]
            if address not in blood_vessel_addresses:
                slow_addresses.append(address)
            else:
                # TODO - avoid conflict
                pass





    return blood_vessel_addresses, None, None, None




def main():

    # LOAD CONFIG
    config = ConfigParser.RawConfigParser()
    if not config.read('config.properties'):
        raise IOError("Config file (config.properties) not found")

    # LOAD PARAMETERS
    parameters = dict()
    # Get all options in parameters section and add to the dictionary
    for i in config.options("ParametersSection"):
        if i == 'max_depth':
            parameters[i] = config.getint("ParametersSection", i)
        else:
            parameters[i] = config.getfloat("ParametersSection", i)

    # LOAD GRID ATTRIBUTES
    total_size = [int(a) for a in config.get("GridSection", "total_shape").split(",")]
    tile_arrangement = [int(a) for a in config.get("GridSection", "tile_arrangement").split(",")]

    # LOAD CELL ATTRIBUTES
    attributes = config.get("CellSection", "attributes").split(",")

    # LOAD RUN PARAMETERS
    time_limit = config.getint("RunParametersSection", "time_limit")
    method = config.get("RunParametersSection", "method")

    # TODO - LOAD INITIALISATION
    blood_vessels = []
    fast_bacteria = []
    slow_bacteria = []
    macrophages = []

    blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise(config)

    topology = TB_Model.TwoDimensionalTopology(tile_arrangement, total_size, attributes, parameters, [[3,3]], [[1,1]], [[9,9]], [[7,1]])

    if method == 'single':
        run_single(topology, time_limit)
    elif method == 'many_serial':
        run_many_serial(topology, time_limit)
    elif method == 'many_parallel':
        run_many_parallel(topology, time_limit)
    else:
        raise Exception("Invalid method")


if __name__ == '__main__':
    main()