import TB_Model

def run_single():
    pass

def run_many_serial():
    pass

def run_many_parallel():
    pass

def main():
    params = dict()
    params['max_depth'] = 3
    params['initial_oxyegn'] = 1.0
    atts = ['oxygen', 'blood_vessel', 'contents']
    topology = TB_Model.TwoDimensionalTopology([2, 2], [10, 10], atts, params, [[3,3]], [[1,1]], [[9,9]], [[7,1]])

    for a in topology.automata:
        print a.grid


if __name__ == '__main__':
    main()