# Python_TB_Model

Tuberculosis Automaton Model
Michael J. Pitcher
School of Computer Science, University of St. Andrews

Creates a parallelised cellular automaton / agent hybrid model, where the grid is split into multiple smaller grids
to improve performance. A series of automata run updates on the tiles and agents which are on the grids, and these
create a series of potential events. Once conflicting events are resolved (in runner.py) the acceptable events
are passed back to automata which update the tiles, and the process repeats as needed.

Some terminology:

Cell - the individual blocks of the tile (here represented as a dictionary of attributes)
Tile - the smaller grids which constitute the larger overall grid
Size - Number of cells in a tile
Shape - Arrangement of cells in a tile
Address - an (x,y,z,etc.) collection of co-ordinates for a cell
Halo - the cells required by a tile for updating which are not a part of the tile (belong to other tiles)
Danger Zone - the cells in a tile which will be required by other tiles (i.e. are part of other tiles' halos)
