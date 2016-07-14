import cProfile
import TB_Model
import numpy as np
from ipyparallel import Client


def setup_topology():

    params = dict()
    params['max_depth'] = 3
    params['caseum_threshold_to_reduce_diffusion'] = 2
    params['oxygen_diffusion_caseum_reduction'] = 1.5
    params['chemotherapy_diffusion_caseum_reduction'] = 1.5
    params['initial_oxygen'] = 1.5
    params['oxygen_diffusion'] = 0.0
    params['chemotherapy_diffusion'] = 0.0
    params['caseum_distance_to_reduce_diffusion'] = 2
    params['spatial_step'] = 0.2
    params['chemotherapy_schedule1_start'] = 99
    params['chemotherapy_schedule2_start'] = 200
    params['oxygen_from_source'] = 0.0
    params['chemokine_diffusion'] = 0.0
    params['chemokine_decay'] = 0.0
    params['chemokine_from_macrophage'] = 0
    params['bacteria_threshold_for_t_cells'] = 100
    params['chemotherapy_scale_for_kill_macrophage'] = 101
    params['oxygen_uptake_from_bacteria'] = 0
    params['chemokine_from_bacteria'] = 0
    params['bacteria_replication_fast_upper'] = 999
    params['bacteria_replication_fast_lower'] = 998
    params['bacteria_replication_slow_upper'] = 999
    params['bacteria_replication_slow_lower'] = 998
    params['chemotherapy_scale_for_kill_fast_bacteria'] = 101
    params['chemotherapy_scale_for_kill_slow_bacteria'] = 101
    params['resting_macrophage_age_limit'] = 999
    params['active_macrophage_age_limit'] = 999
    params['infected_macrophage_age_limit'] = 999
    params['chronically_infected_macrophage_age_limit'] = 999
    params['resting_macrophage_movement_time'] = 1000
    params['active_macrophage_movement_time'] = 1000
    params['infected_macrophage_movement_time'] = 1000
    params['chronically_infected_macrophage_movement_time'] = 1000
    params['prob_resting_macrophage_random_move'] = 0
    params['minimum_chemokine_for_resting_macrophage_movement'] = 0
    params['bacteria_to_turn_chronically_infected'] = 10
    params['chemokine_scale_for_macrophage_activation'] = 0
    params['chemokine_scale_for_macrophage_deactivation'] = 100
    params['time_step'] = 4
    params['oxygen_scale_for_metabolism_change_to_slow'] = 100
    params['oxygen_scale_for_metabolism_change_to_fast'] = 0
    params['macrophage_recruitment_probability'] = 0
    params['chemotherapy_from_source'] = 1
    params['chemotherapy_decay'] = 0.2

    params['bacteria_to_burst_macrophage'] = 20

    atts = ['blood_vessel', 'contents', 'oxygen', 'oxygen_diffusion_rate', 'chemotherapy_diffusion_rate',
        'chemotherapy', 'chemokine']

    blood_vessels = [[1, 1]]
    fast_bacteria = []
    slow_bacteria = []
    macrophages = [[2, 2]]
    topology = TB_Model.TwoDimensionalTopology([2, 2], [100, 100], atts, params, blood_vessels, fast_bacteria,
                                                    slow_bacteria, macrophages)

    return topology


def sort_out_halos(topology):
    dz = []
    for i in topology.automata:
        dz.append(i.get_danger_zone())
    halos = topology.create_halos(dz)
    for i in range(4):
        topology.automata[i].set_halo(halos[i])

    return topology


def run_many_serial(topology, time_limit):

    print "Method: Many serial"
    number_tiles = len(topology.automata)

    # HALOS
    dz_addresses = []
    for tile_id_of_event in range(number_tiles):
        dz_addresses.append(topology.automata[tile_id_of_event].danger_zone_addresses)

    halo_addresses = topology.external_addresses_required

    danger_zone_values = []
    received_events = dict()
    events_to_return = dict()
    addresses_processed = dict()

    for a in range(number_tiles):
        received_events[a] = []
        events_to_return[a] = []
        addresses_processed[a] = []

    for t in range(1, time_limit+1):
        print "TIME-STEP:", t

        max_oxygen = 0.0
        max_chemotherapy = 0.0
        max_chemokine = 0.0
        number_bacteria = 0

        # 1. Get values from engines
        for tile_id in range(number_tiles):
            automaton = topology.automata[tile_id]

            # ---------------- BACK -----------------------------
            danger_zone_values.append(automaton.get_danger_zone())
            max_oxygen = max(max_oxygen, automaton.max_oxygen_local)
            max_chemotherapy = max(max_chemotherapy, automaton.max_chemotherapy_local)
            max_chemokine = max(max_chemokine, automaton.max_chemokine_local)
            number_bacteria += len(automaton.bacteria)
            # ---------------- BACK -----------------------------

        # 2. Construct halos
        halos = construct_halos(topology, dz_addresses, danger_zone_values, halo_addresses)

        # Reset values
        danger_zone_values = []

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

        events_to_return = veto_conflicting_events(received_events, topology)

        # SEND EVENTS TO TILES TO PERFORM
        for tile_id_of_event in range(len(topology.automata)):
            automaton = topology.automata[tile_id_of_event]
            # ------------------- OUT ----------------------------
            # 7 & 8. Send events, perform events
            automaton.process_events(events_to_return[tile_id_of_event])
            # ------------------- OUT ----------------------------
            # print automaton.grid


def run_many_parallel(topology, time_limit, json_controller_path):

    print "Method: Many parallel"
    number_tiles = len(topology.automata)

    # HALOS
    dz_addresses = []
    for tile_id_of_event in range(number_tiles):
        dz_addresses.append(topology.automata[tile_id_of_event].danger_zone_addresses)

    halo_addresses = topology.external_addresses_required

    received_events = dict()
    events_to_return = dict()
    addresses_processed = dict()

    for a in range(number_tiles):
        received_events[a] = []
        events_to_return[a] = []
        addresses_processed[a] = []

    # SET UP CONNECTION TO ENGINES
    direct_view = create_view(json_controller_path)
    print len(direct_view), "engines are running"
    assert len(direct_view) == number_tiles, \
        "Incorrect set-up - number of engines {0} does not match expected number of tiles {1}"\
            .format(len(direct_view), number_tiles)

    import TB_Model
    send_module(direct_view, TB_Model)
    with direct_view.sync_imports():
        import csv
        import numpy

    direct_view.scatter('automaton', topology.automata)
    direct_view.scatter('danger_zone', [[], ] * number_tiles)
    direct_view.scatter('max_oxygen_local', [[], ] * number_tiles)
    direct_view.scatter('max_chemotherapy_local', [[], ] * number_tiles)
    direct_view.scatter('max_chemokine_local', [[], ] * number_tiles)
    direct_view.scatter('number_bacteria_local', [[], ] * number_tiles)
    direct_view.scatter('potential_events', [[], ] * number_tiles)

    for t in range(1, time_limit + 1):
        print "TIME-STEP:", t

        # 1. Get values from engines
        # ---------------- BACK -----------------------------
        direct_view.apply(parallel_update_engine_variables())

        danger_zone_values = direct_view.gather('danger_zone')
        local_oxygen_maxima = direct_view.gather('max_oxygen_local')
        max_oxygen = max(local_oxygen_maxima)
        local_chemotherapy_maxima = direct_view.gather('max_chemotherapy_local')
        max_chemotherapy = max(local_chemotherapy_maxima)
        local_chemokine_maxima = direct_view.gather('max_chemokine_local')
        max_chemokine = max(local_chemokine_maxima)
        bacteria_counts = direct_view.gather('number_bacteria_local')
        number_bacteria = sum(bacteria_counts)
        # ---------------- BACK -----------------------------

        # 2. Construct halos
        halos = construct_halos(topology, dz_addresses, danger_zone_values, halo_addresses)

        # 3. Send halos & values out to engines
        # ---------------- OUT -----------------------------
        direct_view.scatter('halo', halos)
        direct_view.push(dict(max_oxygen_global=max_oxygen, max_chemotherapy_global=max_chemotherapy,
                              max_chemokine_global=max_chemokine, number_bacteria_global=number_bacteria))
        # 4. Calculate cell values
        direct_view.apply(parallel_run_updates())
        # ---------------- OUT -----------------------------

        # ----------------- BACK ----------------------------
        events = direct_view.gather('potential_events')

        total_num_events = 0
        automata_with_events_left = []

        for tile_id in range(len(events)):
            if len(events[tile_id]) > 0:
                received_events[tile_id] = events[tile_id]
                total_num_events += len(events[tile_id])
                automata_with_events_left.append(tile_id)
            else:
                received_events[tile_id] = []
        # ----------------- BACK ----------------------------

        events_to_return = veto_conflicting_events(received_events, topology)

        # 7. SEND EVENTS TO TILES TO PERFORM
        direct_view.scatter('accepted_events', events_to_return.values())
        # 8. PERFORM EVENTS
        direct_view.apply(parallel_process_events())

    automata = direct_view.gather('automaton')


def create_view(path):
    """Create a (blocking) direct view to the processing engines"""
    # Create a connection to the client
    if path == "None":
        client = Client()
    else:
        client = Client(path)
    # Create direct view
    direct_view = client[:]
    # Blocking on
    direct_view.block = True
    return direct_view


def veto_conflicting_events(received_events, topology):

    total_num_events = 0
    automata_with_events_left = []
    addresses_processed = dict()
    events_to_return = dict()

    for tile_id in range(len(received_events)):
        addresses_processed[tile_id] = []
        events_to_return[tile_id] = []
        if len(received_events[tile_id]) > 0:
            automata_with_events_left.append(tile_id)
            total_num_events += len(received_events[tile_id])

    for i in range(total_num_events):
        # Pick a tile at random to pull events from
        # TODO - COMP - check validity of this approach
        tile_id_of_event = automata_with_events_left[np.random.randint(0, len(automata_with_events_left))]

        event = received_events[tile_id_of_event].pop()

        # Wholly internal
        if event.internal:
            flag = True
            for address in event.dependant_addresses:
                if address in addresses_processed[tile_id_of_event]:
                    flag = False
                    break
            if flag:
                # Loop through the impacted addresses. If an impacted address has already been processed, then
                # it will not be affected
                for impacted_address in event.impacted_addresses_potential:
                    if impacted_address not in addresses_processed[tile_id_of_event]:
                        event.impacted_addresses_allowed.append(impacted_address)
                        addresses_processed[tile_id_of_event].append(impacted_address)
                events_to_return[tile_id_of_event].append(event)
        else:  # Crosses a boundary
            flag = True
            for address in event.dependant_addresses:
                global_address = topology.local_to_global(tile_id_of_event, address)
                local_tile_id, local_address = topology.global_to_local(global_address)
                if local_address in addresses_processed[local_tile_id]:
                    flag = False
                    break
            # Event is acceptable
            if flag:
                # Determine which of the impacted addresses can be processed
                tiles_impacted = []
                for impacted_address in event.impacted_addresses_potential:
                    global_address = topology.local_to_global(tile_id_of_event, impacted_address)
                    local_tile_id, local_address = topology.global_to_local(global_address)
                    if local_address not in addresses_processed[local_tile_id]:
                        event.impacted_addresses_allowed.append(impacted_address)
                        addresses_processed[local_tile_id].append(local_address)
                        tiles_impacted.append(local_tile_id)
                # Having determine which addresses can be impacted and what tiles are affected, sort out events
                for impacted_tile_id in tiles_impacted:
                    # If it's the original tile, just return the event and record addresses
                    if impacted_tile_id == tile_id_of_event:
                        events_to_return[tile_id_of_event].append(event)
                        addresses_processed[tile_id_of_event] += event.impacted_addresses
                    # If this event impacts a new tile, need to clone the event and add addresses relative to
                    # the new tile
                    else:
                        new_addresses = []
                        for address in event.impacted_addresses_allowed:
                            new_addresses.append(
                                topology.local_to_local(tile_id_of_event, address, impacted_tile_id))
                        new_event = event.clone(new_addresses)
                        new_event.impacted_addresses_allowed = new_addresses
                        events_to_return[impacted_tile_id].append(new_event)

        if len(received_events[tile_id_of_event]) == 0:
            automata_with_events_left.remove(tile_id_of_event)

    return events_to_return


def parallel_update_engine_variables():

    def function():
        max_oxygen_local[0] = automaton[0].max_oxygen_local
        max_chemotherapy_local[0] = automaton[0].max_chemotherapy_local
        max_chemokine_local[0] = automaton[0].max_chemokine_local

        danger_zone[0] = automaton[0].get_danger_zone()

        number_bacteria_local[0] = len(automaton[0].bacteria)

    return function


def parallel_run_updates():

    def function():
        automaton[0].set_halo(halo[0])
        automaton[0].set_max_oxygen_global(max_oxygen_global)
        automaton[0].set_max_chemotherapy_global(max_chemotherapy_global)
        automaton[0].set_max_chemokine_global(max_chemokine_global)
        automaton[0].set_global_bacteria_number(number_bacteria_global)

        automaton[0].update()
        potential_events[0] = automaton[0].potential_events

    return function


def parallel_process_events():

    def function():
        automaton[0].process_events(accepted_events[0])

    return function


def send_module(view, mod):
    """Send the module to the engines for import

    Given a module, and a ipyparallel direct view, copies the pyc file of module onto the
    working directory of the engines, so the module can be imported
    """
    # TODO - this isn't very good and there *must* be a better way to pass modules to the engines

    # Open the pyc file of the class and read the data
    file_path = mod.__file__
    with open(file_path) as f:
        data = f.read()

    # Function to be run on each remote engine
    def _write_module(filename, file_data):

        # Get a new file location - current working directory + filename
        import os
        # Split as os.path is absolute and only want the filename
        filename = os.getcwd() + '/' + os.path.split(filename)[1]

        # Create a new file using the filename and write the original data to it
        with open(filename, 'w') as file_:
            file_.write(file_data)
        return filename

    # Run the function on every engine on the (slight) off-chance they have different working directories
    return view.apply(_write_module, file_path, data)


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
    topology = setup_topology()
    sort_out_halos(topology)

    pr = cProfile.Profile()

    pr.enable()
    run_many_serial(topology, 10)
    pr.disable()

    pr.print_stats(sort='tottime')


if __name__ == '__main__':
    main()