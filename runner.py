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


def main():

    config = ConfigParser.RawConfigParser()
    if not config.read('config.properties'):
        raise IOError("Config file (config.properties) not found")

    parameters = dict()
    # Get all options in parameters section and add to the dictionary
    for i in config.options("ParametersSection"):
        if i == 'max_depth':
            parameters[i] = config.getint("ParametersSection", i)
        else:
            parameters[i] = config.getfloat("ParametersSection", i)

    attributes = config.get("CellSection", "attributes").split(",")

    time_limit = config.getint("RunParametersSection", "time_limit")
    method = config.get("RunParametersSection", "method")

    # TODO - initialise

    topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], attributes, parameters, [[3,3]], [[1,1]], [[9,9]], [[7,1]])

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