import ConfigParser
import itertools
import numpy as np
import time
from ipyparallel import Client

#import pyximport; pyximport.install(pyimport=True)
import TB_Model


def run_single(topology, time_limit):
    pass

    print "Method: Single"

    halo_addresses = topology.external_addresses_required

    empty_halo = dict()
    for a in halo_addresses:
        empty_halo[a] = None

    print 'Running simulation...'
    simulation_start_time = time.time()
    for t in range(1, time_limit + 1):
        topology.automata[0].set_halo(empty_halo)
        print "TIME-STEP:", t
        automaton = topology.automata[0]
        automaton.set_max_oxygen_global(automaton.max_oxygen_local)
        automaton.set_max_chemotherapy_global(automaton.max_chemotherapy_local)
        automaton.set_max_chemokine_global(automaton.max_chemokine_local)
        automaton.set_global_bacteria_number(len(automaton.bacteria))

        automaton.update()

        events = [automaton.potential_events]

        events_to_return = veto_conflicting_events(events, topology)

        automaton.process_events(events_to_return[0])
    simulation_end_time = time.time()
    print 'Simulation complete. Time taken:', simulation_end_time-simulation_start_time


def run_many_serial(topology, time_limit):

    print "Method: Many serial"

    number_tiles = len(topology.automata)

    # HALOS
    dz_addresses = []
    for tile_id_of_event in range(number_tiles):
        dz_addresses.append(topology.automata[tile_id_of_event].danger_zone_addresses)

    danger_zone_values = []
    received_events = dict()
    events_to_return = dict()
    addresses_processed = dict()

    for a in range(number_tiles):
        received_events[a] = []
        events_to_return[a] = []
        addresses_processed[a] = []

    print 'Running simulation...'
    simulation_start_time = time.time()
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
        halos = construct_halos(topology, danger_zone_values)

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

    simulation_end_time = time.time()
    print 'Simulation complete. Time taken:', simulation_end_time - simulation_start_time


def run_many_parallel(topology, time_limit, json_controller_path):

    print "Method: Many parallel"
    number_tiles = len(topology.automata)

    # HALOS
    dz_addresses = []
    for tile_id_of_event in range(number_tiles):
        dz_addresses.append(topology.automata[tile_id_of_event].danger_zone_addresses)

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
    print 'Send code to engines'
    send_module(direct_view, TB_Model)
    with direct_view.sync_imports():
        import csv
        import numpy

    print 'Sending automata to engines'
    send_start_time = time.time()
    direct_view.scatter('automaton', topology.automata)
    send_end_time = time.time()
    print 'Automata sent. Time taken:', send_end_time-send_start_time

    # Scatter empty lists
    direct_view.scatter('stats', [dict(), ] * number_tiles)
    direct_view.scatter('potential_events', [[], ] * number_tiles)
    direct_view.apply(parallel_update_engine_variables())

    print 'Running simulation'
    simulation_start_time = time.time()
    for t in range(1, time_limit + 1):
        print "TIME-STEP:", t

        # 1. Get values from engines
        # ---------------- BACK -----------------------------
        stats = direct_view.gather('stats')
        danger_zone_values = []
        max_oxygen = 0.0
        max_chemotherapy = 0.0
        max_chemokine = 0.0
        number_bacteria = 0
        for i in stats:
            danger_zone_values.append(i['danger_zone'])
            max_oxygen = max(max_oxygen, i['max_oxygen_local'])
            max_chemotherapy = max(max_chemotherapy, i['max_chemotherapy_local'])
            max_chemokine = max(max_chemokine, i['max_chemokine_local'])
            number_bacteria += i['number_bacteria_local']
        # ---------------- BACK -----------------------------

        # 2. Construct halos
        halos = construct_halos(topology, danger_zone_values)

        # 3. Send halos & values out to engines
        # ---------------- OUT -----------------------------
        direct_view.scatter('halo', halos)
        # 4. Calculate cell values
        direct_view.apply(parallel_run_updates(), max_oxygen, max_chemotherapy, max_chemokine, number_bacteria)
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

    simulation_end_time = time.time()
    print 'Simulation complete. Time taken:', simulation_end_time-simulation_start_time


def construct_halos(topology, danger_zone_values):

    return topology.create_halos(danger_zone_values)


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


def initialise(config, total_shape):

    available_addresses = []

    # TODO - only works for 2d
    for a in itertools.product(range(total_shape[0]), range(total_shape[1])):
        available_addresses.append(a)

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")
    blood_vessel_addresses = []

    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", "blood_vessels_hard_code").split("/")
        for b in bv_list:
            address = tuple(int(c) for c in b.split(","))
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
            address = tuple(int(c) for c in a.split(","))
            if address in available_addresses:
                fast_addresses.append(address)
                available_addresses.remove(address)
            else:
                # TODO - avoid conflict
                pass

        slow_list = config.get("InitialiseSection", "bacteria_slow_hard_code").split("/")
        assert len(available_addresses) > len(slow_list)
        for a in slow_list:
            address = tuple(int(c) for c in a.split(","))
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
                        addresses_processed[tile_id_of_event] += event.impacted_addresses_allowed
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
        stats[0]['danger_zone'] = automaton[0].get_danger_zone()
        stats[0]['max_oxygen_local'] = automaton[0].max_oxygen_local
        stats[0]['max_chemotherapy_local'] = automaton[0].max_chemotherapy_local
        stats[0]['max_chemokine_local'] = automaton[0].max_chemokine_local
        stats[0]['number_bacteria_local'] = len(automaton[0].bacteria)

    return function


def parallel_run_updates():

    def function(max_oxygen_global, max_chemotherapy_global, max_chemokine_global, number_bacteria_global):
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
        stats[0]['danger_zone'] = automaton[0].get_danger_zone()
        stats[0]['max_oxygen_local'] = automaton[0].max_oxygen_local
        stats[0]['max_chemotherapy_local'] = automaton[0].max_chemotherapy_local
        stats[0]['max_chemokine_local'] = automaton[0].max_chemokine_local
        stats[0]['number_bacteria_local'] = len(automaton[0].bacteria)

    return function


def main():

    print '------------------------'
    print 'TB Simulation Automaton'
    print '------------------------'

    whole_start_time = time.time()

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

    # LOAD INITIALISATION
    blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise(config, total_shape)

    print 'Constructing topology...'
    construction_start_time = time.time()
    topology = TB_Model.TwoDimensionalTopology(tile_arrangement, total_shape, attributes, parameters, blood_vessels,
                                               fast_bacteria, slow_bacteria, macrophages)
    construction_end_time = time.time()
    print 'Complete. Time taken for construction: ', construction_end_time-construction_start_time

    if method == 'single':
        run_single(topology, time_limit)
    elif method == 'many_serial':
        run_many_serial(topology, time_limit)
    elif method == 'many_parallel':
        # Remote section
        json_controller_path = config.get("RemoteSection", "json_controller_path")
        run_many_parallel(topology, time_limit, json_controller_path)
    else:
        raise Exception("Invalid method")

    whole_end_time = time.time()

    print 'Whole process time:', whole_end_time-whole_start_time

if __name__ == '__main__':
    main()