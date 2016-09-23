import numpy as np
import TB_Model_v2
import ConfigParser

def initialise(config, total_shape):

    available_addresses = []

    # TODO - only works for 2d
    for a in itertools.product(range(total_shape[0]), range(total_shape[1])):
        available_addresses.append(a)

    # BLOOD_VESSELS
    blood_vessels_method = config.get("InitialiseSection", "blood_vessels")
    blood_vessel_addresses = []

    if blood_vessels_method == 'hard_code':
        bv_list = config.get("InitialiseSection", 'blood_vessels_hard_code').split('/')
        for b in bv_list:
            address = tuple(int(c) for c in b.split(","))
            available_addresses.remove(address)
            blood_vessel_addresses.append(address)
    elif blood_vessels_method == 'from_file':
        path = config.get("InitialiseSection", 'blood_vessels_from_file')
        # Add the values to the list
        list_of_vessels = [line.rstrip('\n') for line in open(path)]
        for index in range(len(list_of_vessels)):
            if float(list_of_vessels[index]) > 0.0:
                try:
                    address = np.unravel_index(index, total_shape)
                except ValueError:
                    raise Exception(index, total_shape)
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


def main():
    print '------------------------'
    print 'TB Simulation Automaton'
    print '------------------------'

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

    # LOAD CELL ATTRIBUTES
    attributes = config.get("CellSection", "attributes").split(",")

    # LOAD RUN PARAMETERS
    time_limit = int(config.getint("RunParametersSection", "time_limit") / parameters['time_step'])

    # LOAD INITIALISATION
    blood_vessels, fast_bacteria, slow_bacteria, macrophages = initialise(config, total_shape)

    print 'Constructing topology...'

    automaton = TB_Model_v2.Automaton(total_shape, parameters, blood_vessels, fast_bacteria, slow_bacteria, macrophages)

    automaton.run()

    print 'Completed'

if __name__ == '__main__':
    main()