import TB_Model
import ConfigParser
import itertools
import numpy as np


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

    received_events = dict()
    events_to_return = dict()
    addresses_processed = dict()

    for a in range(number_tiles):
        received_events[a] = []
        events_to_return[a] = []
        addresses_processed[a] = []

    for t in range(time_limit):
        print "TIME-STEP:", t

        max_oxygen = 0.0
        max_chemotherapy = 0.0
        max_chemokine = 0.0
        number_bacteria = 0

        # 1. Get DZ values from engines
        for tile_id in range(number_tiles):
            automaton = topology.automata[tile_id]

            # ---------------- BACK -----------------------------
            values.append(automaton.get_danger_zone())
            max_oxygen = max(max_oxygen, automaton.max_oxygen_local)
            max_chemotherapy = max(max_chemotherapy, automaton.max_chemotherapy_local)
            max_chemokine = max(max_chemokine, automaton.max_chemokine_local)
            number_bacteria += len(automaton.bacteria)
            # ---------------- BACK -----------------------------

        # 2. Construct halos
        halos = construct_halos(topology, dz_addresses, values, halo_addresses)

        # Reset values
        values = []

        automata_with_events_left = []
        total_num_events = 0
        for automaton in topology.automata:
            # Reset
            events_to_return[automaton.tile_id] = []
            addresses_processed[automaton.tile_id] = []

            # UPDATE
            # ------------------ OUT ----------------------------
            # 3. Send halos to engines
            automaton.set_halo(halos[automaton.tile_id])
            automaton.set_max_oxygen_global(max_oxygen)
            automaton.set_max_chemotherapy_global(max_chemotherapy)
            automaton.set_max_chemokine_global(max_chemokine)
            automaton.set_global_bacteria_number(number_bacteria)
            # 4. Calculate cell values
            automaton.update()
            # ------------------ OUT ----------------------------

            # ----------------- BACK ----------------------------
            events = automaton.potential_events

            if len(events) > 0:
                # 5. Engines return list of potential events
                received_events[automaton.tile_id] = events
                total_num_events += len(events)
                automata_with_events_left.append(automaton.tile_id)
            else:
                received_events[automaton.tile_id] = []
            # ----------------- BACK ----------------------------

        # 6. VETO CONFLICTING EVENTS
        for i in range(total_num_events):
            tile_id = automata_with_events_left[np.random.randint(0, len(automata_with_events_left))]

            event = received_events[tile_id].pop()

            # Wholly internal
            if event.internal:
                flag = True
                for address in event.addresses_affected:
                    if address in addresses_processed[tile_id]:
                        flag = False
                        break
                if flag:
                    events_to_return[tile_id].append(event)
                    addresses_processed[tile_id] += event.addresses_affected
            else:  # Crosses a boundary
                flag = True
                tiles_affected = []
                for address in event.addresses_affected:
                    global_address = topology.local_to_global(tile_id, address)
                    local_tile_id, local_address = topology.global_to_local(global_address)
                    tiles_affected.append(local_tile_id)
                    if local_address in addresses_processed[local_tile_id]:
                        flag = False
                        break
                if flag:
                    for index in range(len(tiles_affected)):

                        tile_affected = tiles_affected[index]

                        if tile_affected == tile_id:
                            events_to_return[tile_id].append(event)
                            addresses_processed[tile_id] += event.addresses_affected
                        else:
                            new_addresses = []
                            for address in event.addresses_affected:
                                new_addresses.append(topology.local_to_local(tile_id, address, tile_affected))
                            new_event = event.clone(new_addresses)
                            events_to_return[tile_affected].append(new_event)
                            addresses_processed[tile_affected] += new_addresses

            if len(received_events[tile_id]) == 0:
                automata_with_events_left.remove(tile_id)

        # SEND EVENTS TO TILES TO PERFORM
        for tile_id in range(len(topology.automata)):
            automaton = topology.automata[tile_id]
            # ------------------- OUT ----------------------------
            # 7 & 8. Send events, perform events
            automaton.process_events(events_to_return[tile_id])
            # ------------------- OUT ----------------------------
            # print automaton.grid


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


def initialise(config, total_shape):

    available_addresses = []

    # TODO - only works for 2d
    for a in itertools.product(range(total_shape[0]), range(total_shape[1])):
        available_addresses.append(list(a))

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")
    blood_vessel_addresses = []

    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", "blood_vessels_hard_code").split("/")
        for b in bv_list:
            address = [int(c) for c in b.split(",")]
            available_addresses.remove(address)
            blood_vessel_addresses.append(address)
    elif blood_vessels_method == 'random':
        number = config.getint("InitialiseSection", "blood_vessels_random_number")
        assert len(available_addresses) > number
        for i in range(number):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            blood_vessel_addresses.append(address)

    # FAST BACTERIA
    bacteria_method = config.get("InitialiseSection","bacteria")
    fast_addresses = []
    slow_addresses = []

    if bacteria_method == 'hard_code':
        fast_list = config.get("InitialiseSection", "bacteria_fast_hard_code").split("/")
        assert len(available_addresses) > len(fast_list)
        for a in fast_list:
            address = [int(c) for c in a.split(",")]
            if address in available_addresses:
                fast_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass

        slow_list = config.get("InitialiseSection", "bacteria_slow_hard_code").split("/")
        assert len(available_addresses) > len(slow_list)
        for a in slow_list:
            address = [int(c) for c in a.split(",")]
            if address not in blood_vessel_addresses:
                slow_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass
    elif bacteria_method == 'random':
        number_fast = config.getint("InitialiseSection", "bacteria_fast_random_number")
        assert len(available_addresses) > number_fast
        for i in range(number_fast):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            fast_addresses.append(address)

        number_slow = config.getint("InitialiseSection", "bacteria_slow_random_number")
        assert len(available_addresses) > number_slow
        for i in range(number_slow):
            address = available_addresses.pop(np.random.randint(0, len(available_addresses)))
            slow_addresses.append(address)

    # MACROPHAGES
    macrophage_method = config.get("InitialiseSection", "macrophages")
    macrophage_addresses = []
    if macrophage_method == 'random':
        number = config.getint("InitialiseSection", "macrophages_random_number")

        # Make sure there's enough room
        assert len(available_addresses) > number

        for i in range(number):
            address = available_addresses.pop(np.random.randint(0,len(available_addresses)))
            macrophage_addresses.append(address)

    return blood_vessel_addresses, fast_addresses, slow_addresses, macrophage_addresses


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
    total_shape = [int(a) for a in config.get("GridSection", "total_shape").split(",")]
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

    blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise(config, total_shape)

    topology = TB_Model.TwoDimensionalTopology(tile_arrangement, total_shape, attributes, parameters, blood_vessels,
                                               fast_bacteria, slow_bacteria, macrophages)

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