import TB_Model


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

    # TODO - config

    params = dict()
    params['max_depth'] = 3
    params['initial_oxygen'] = 1.0
    atts = ['oxygen', 'blood_vessel', 'contents']
    topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, [[3,3]], [[1,1]], [[9,9]], [[7,1]])

    time_limit = 10

    run_many_serial(topology, time_limit)


if __name__ == '__main__':
    main()