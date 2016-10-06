import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import gridspec


class Displayer:

    def __init__(self, output_directory, vessels, grid_shape):
        self.output_directory = output_directory
        self.vessels = vessels
        self.shape = grid_shape

    def display_numbers(self):
        pass

    def grid_animation(self, movie_name, interval=400, legend=False, display=True):

        print "Collecting data..."
        size = reduce(lambda i, j: i * j, self.shape)
        with open(self.output_directory + "/data_test.txt") as f:
            cell_data = [float(x.strip('\n')) for x in f.readlines()]
        with open(self.output_directory + "/oxygen_test.txt") as f:
            oxygen_data = [float(x.strip('\n')) for x in f.readlines()]
        with open(self.output_directory + "/chemo1.txt") as f:
            chemotherapy_data = [float(x.strip('\n')) for x in f.readlines()]
        with open(self.output_directory + "/ckine.txt") as f:
            chemokine_data = [float(x.strip('\n')) for x in f.readlines()]

        assert len(cell_data) % size == len(oxygen_data) % size == len(chemotherapy_data) % size \
            == len(chemokine_data) % size == 0.0

        time_steps = len(oxygen_data) / size

        # GATHER DATA
        contents_grids = []
        for t in range(time_steps):
            contents_grids.append(np.array(cell_data[0 + (t * size):size + (t * size)]).reshape(self.shape))
        fast_bacs = []
        fast_rest_bacs = []
        slow_bacs = []
        slow_rest_bacs = []
        rest_macs = []
        active_macs = []
        inf_macs = []
        chr_inf_macs = []
        t_cells = []
        caseum = []

        for time_step in range(len(contents_grids)):
            fast_bacs.append([])
            fast_rest_bacs.append([])
            slow_bacs.append([])
            slow_rest_bacs.append([])
            rest_macs.append([])
            active_macs.append([])
            inf_macs.append([])
            chr_inf_macs.append([])
            t_cells.append([])
            caseum.append([])

            grid = contents_grids[time_step]
            for y in range(grid.shape[0]):
                for x in range(grid.shape[1]):
                    # FAST BAC
                    if grid[(x, y)] == 1.0:
                        fast_bacs[time_step].append((x, y))
                    # FAST BAC REST
                    elif grid[(x, y)] == 1.25:
                        fast_rest_bacs[time_step].append((x, y))
                    # SLOW BAC
                    elif grid[(x, y)] == 2.0:
                        slow_bacs[time_step].append((x, y))
                    # SLOW BAC REST
                    elif grid[(x, y)] == 2.25:
                        slow_rest_bacs[time_step].append((x, y))
                    # REST MAC
                    elif grid[(x, y)] == 4.0:
                        rest_macs[time_step].append((x, y))
                    # ACTIVE MAC
                    elif grid[(x, y)] == 5.0:
                        active_macs[time_step].append((x, y))
                    # INF MAC
                    elif grid[(x, y)] == 6.0:
                        inf_macs[time_step].append((x, y))
                    # CHR INF MAC
                    elif grid[(x, y)] == 7.0:
                        chr_inf_macs[time_step].append((x, y))
                    # T CELL
                    elif grid[(x, y)] == 3.0:
                        t_cells[time_step].append((x, y))
                    # CASEUM
                    elif grid[(x, y)] == 100.0:
                        caseum[time_step].append((x, y))

        def update_plot(time_step):
            plt.clf()
            if legend:
                gs = gridspec.GridSpec(2, 1, height_ratios=[6, 1])
                plt.subplot(gs[0])

            plt.axis([0, self.shape[0], self.shape[1], 0])
            plt.xticks([])
            plt.yticks([])
            plt.suptitle("TB Automaton", fontsize=14, fontweight='bold')
            plt.title('Time = ' + str(time_step) + " hours", fontsize=10)

            bv = plt.scatter([v[1] for v in self.vessels], [v[0] for v in self.vessels],
                        s=20, color='red', marker="D")  # RED
            fb = plt.scatter([fb[1] for fb in fast_bacs[time_step]], [fb[0] for fb in fast_bacs[time_step]],
                        s=1, color='#0F63AE')  # BLUE
            frb = plt.scatter([fbr[1] for fbr in fast_rest_bacs[time_step]], [fbr[0] for fbr in fast_rest_bacs[time_step]],
                        s=1, color='#0A4579')  # DEEP BLUE
            sb = plt.scatter([sb[1] for sb in slow_bacs[time_step]], [sb[0] for sb in slow_bacs[time_step]],
                        s=1, color='#851f98', marker="D")  # PURPLE
            srb = plt.scatter([sbr[1] for sbr in slow_rest_bacs[time_step]], [sbr[0] for sbr in slow_rest_bacs[time_step]],
                        s=1, color='#490746', marker="D")  # DEEP PURPLE
            rm = plt.scatter([rm[1] for rm in rest_macs[time_step]], [rm[0] for rm in rest_macs[time_step]],
                        color='#168964', marker=(5, 1))  # GREEN
            am = plt.scatter([am[1] for am in active_macs[time_step]], [am[0] for am in active_macs[time_step]],
                        color='#00ff45', marker=(5, 1))  # BRIGHT GREEN
            im = plt.scatter([im[1] for im in inf_macs[time_step]], [im[0] for im in inf_macs[time_step]],
                        color='#F1BC41', marker=(5, 1))  # GOLD
            cim = plt.scatter([cim[1] for cim in chr_inf_macs[time_step]],[cim[0] for cim in chr_inf_macs[time_step]],
                        color='#77643a', marker=(5, 1))  # BROWN
            tc = plt.scatter([tc[1] for tc in t_cells[time_step]], [tc[0] for tc in t_cells[time_step]],
                        color='#f9c7ed')  # PINK
            ca = plt.scatter([c[1] for c in caseum[time_step]], [c[0] for c in caseum[time_step]],
                        color='#000000')  # BLACK

            if legend:
                plt.subplot(gs[1])
                plt.xticks([])
                plt.yticks([])
                plt.legend((bv, fb, frb, sb, srb, rm, am, im, cim, tc, ca), ("Blood vessel", "Fast bacterium", "Fast resting bacterium",
                    "Slow bacterium", "Slow resting bacterium", "Resting macrophage", "Active macrophage",
                    "Infected macrophage", "Chr. Infected macrophage", "T-cell", "Caseum"), scatterpoints=1,
                    loc='center', ncol=3, fontsize=8)

        # DISPLAY
        print "Creating animation..."
        fig = plt.figure()
        ani = animation.FuncAnimation(fig, update_plot, frames=xrange(time_steps), interval=interval, blit=False)
        ani.save(output_location + "/" + movie_name + ".mp4", writer='ffmpeg_file')
        if display:
            plt.show()


if __name__ == '__main__':
    # Manual data input
    output_location = '../../../../Comparison/DEBUG/RUTH/1'
    bv_file = '../../../../Comparison/Vessel_files/initialvessel1.txt'
    movie_filename = "TBModel"
    shape = [101, 101]
    with open(bv_file) as bv_file:
        bvs = [float(line.strip('\n')) for line in bv_file.readlines()]
    integer_locations = [il for il in range(len(bvs)) if bvs[il] > 0]
    bv_addresses = [np.unravel_index(a, shape) for a in integer_locations]

    d = Displayer(output_location,bv_addresses,shape)
    d.grid_animation(movie_filename, legend=False)