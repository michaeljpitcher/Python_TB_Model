import TB_Model
import ConfigParser


def run_single():
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

        for tile_id in range(number_tiles):
            # ---------------- BACK -----------------------------
            values.append(topology.automata[tile_id].get_danger_zone())
            # ---------------- BACK -----------------------------

        halos = construct_halos(topology, dz_addresses, values, halo_addresses)

        for automaton in topology.automata:
            print "Running automata: ", automaton.tile_id

            # UPDATE
            # ------------------ OUT ----------------------------
            automaton.set_halo(halos[automaton.tile_id])
            # ------------------ OUT ----------------------------


def run_many_parallel():
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
        if i == 'max_depth' or i == 'time_limit':
            parameters[i] = config.getint("ParametersSection", i)
        else:
            parameters[i] = config.getfloat("ParametersSection", i)

    attributes = config.get("CellSection", "attributes").split(",")

    # TODO - initialise

    topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], attributes, parameters, [[3,3]], [[1,1]], [[9,9]], [[7,1]])

    run_many_serial(topology, parameters['time_limit'])


if __name__ == '__main__':
    main()