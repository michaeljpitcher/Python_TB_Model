import numpy as np
import math

class Tile:

    def __init__(self, shape, attributes):

        self.shape = shape
        self.attributes = attributes
        self.size = reduce(lambda x, y: x * y, shape)

        self.grid = self.create_grid(shape, attributes)

    def create_grid(self, shape, attributes):

        cells = []
        for i in range(self.size):
            cell = dict()
            for att in attributes:
                cell[att] = 0.0
            cells.append(cell)

        grid = np.array(cells).reshape(shape)

        return grid

    def create_work_grid(self):
        self.work_grid = self.grid.copy()
