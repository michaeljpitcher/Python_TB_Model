import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def display_numbers(output_directory):
    pass


def display_grid(output_directory, vessels, shape):

    size = reduce(lambda x, y: x * y, shape)
    with open(output_directory + "/0_data_test.txt") as f:
        cell_data = [float(x.strip('\n')) for x in f.readlines()]
    with open(output_directory + "/0_oxygen_test.txt") as f:
        oxygen_data = [float(x.strip('\n')) for x in f.readlines()]
    with open(output_directory + "/0_chemo1.txt") as f:
        chemotherapy_data = [float(x.strip('\n')) for x in f.readlines()]
    with open(output_directory + "/0_ckine.txt") as f:
        chemokine_data = [float(x.strip('\n')) for x in f.readlines()]

    assert len(oxygen_data) % size == len(chemotherapy_data) % size == len(chemokine_data) % size == 0.0

    timesteps = len(oxygen_data) / size

    contents_grids = []
    for t in range(timesteps):
        contents_grids.append(np.array(cell_data[0 + (t * size):size + (t * size)]).reshape(shape))

    # PRINT
    plt.scatter([v[1] for v in vessels], [v[0] for v in vessels], s=20, color='red', marker="D")

    for y in range(shape[0]):
        for x in range(shape[1]):
            # FAST BAC
            if contents_grids[0][(x,y)] == 1.0:
                plt.scatter([y],[x], s=1, color='#0F63AE')  # BLUE
            # FAST BAC REST
            elif contents_grids[0][(x,y)] == 1.25:
                plt.scatter([y],[x], s=1, color='#0A4579')  # DEEP BLUE
            # SLOW BAC
            elif contents_grids[0][(x,y)] == 2.0:
                plt.scatter([y],[x], s=1, color='#851f98')  # PURPLE
            # SLOW BAC REST
            elif contents_grids[0][(x,y)] == 2.25:
                plt.scatter([y],[x], s=1, color='#490746')  # DEEP PURPLE
            # REST MAC
            elif contents_grids[0][(x,y)] == 4.0:
                plt.scatter([y],[x], color='#168964', marker=(5,1))  # GREEN
            # ACTIVE MAC
            elif contents_grids[0][(x, y)] == 5.0:
                plt.scatter([y],[x], color='#00ff45', marker=(5,1))  # BRIGHT GREEN
            # INF MAC
            elif contents_grids[0][(x, y)] == 6.0:
                plt.scatter([y],[x], color='##eaec17', marker=(5,1))  # YELLOW
            # CHR INF MAC
            elif contents_grids[0][(x, y)] == 7.0:
                plt.scatter([y],[x], color='#c7a30c', marker=(5,1))  # GOLD
            # T CELL
            elif contents_grids[0][(x,y)] == 3.0:
                plt.scatter([y],[x], color='#f9c7ed')  # PINK
            # CASEUM
            elif contents_grids[0][(x, y)] == 100.0:
                plt.scatter([y],[x], color='#000000')  # BLACK

    plt.axis([0, shape[0], shape[1], 0])
    plt.xticks([])
    plt.yticks([])
    plt.title('TB Automaton')
    plt.show()


def main():

    output_location = 'output'
    bv_file = 'Vessel_files/initialvessel1.txt'
    shape = [101,101]

    with open(bv_file) as f:
        bvs = [float(x.strip('\n')) for x in f.readlines()]
    integer_locations = [i for i in range(len(bvs)) if bvs[i] > 0]
    bv_addresses = [np.unravel_index(a, shape) for a in integer_locations]

    display_grid(output_location,bv_addresses,shape)

if __name__ == '__main__':
    main()