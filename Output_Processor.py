import os
import shutil

file_dir = ''
number_of_tiles = 4
tile_shape = (50,50)

if os.path.exists('output'):
    shutil.rmtree('output')
os.makedirs('output')

activemacs = []
caseations = []
chroninfectedmacs = []
infectedmacs = []
intra_bacs = []
restingmacs = []
totals = []
total_cell_tests = []
type1s = []
type1_rs = []
type2s = []
type2_rs = []
type3s = []

contents = []
oxygens = []
chemos = []
ckines = []

for n in range(number_of_tiles):
    # ACTIVEMAC
    activemac_file = [line.rstrip('\n') for line in open(str(n) + '_activemac.txt')]
    activemacs.append(activemac_file)
    # CASEATION
    caseation_file = [line.rstrip('\n') for line in open(str(n) + '_caseation.txt')]
    caseations.append(caseation_file)
    # CHRON_INFECTED_MAC
    chroninfectedmac_file = [line.rstrip('\n') for line in open(str(n) + '_chroninfectedmac.txt')]
    chroninfectedmacs.append(chroninfectedmac_file)
    # INFECTED MAC
    infectedmac_file = [line.rstrip('\n') for line in open(str(n) + '_infectedmac.txt')]
    infectedmacs.append(infectedmac_file)
    # INTRA BAC
    intra_bac_file = [line.rstrip('\n') for line in open(str(n) + '_intra_bac.txt')]
    intra_bacs.append(intra_bac_file)
    # RESTING MACS
    restingmac_file = [line.rstrip('\n') for line in open(str(n) + '_restingmac.txt')]
    restingmacs.append(restingmac_file)
    # TOTAL
    total_file = [line.rstrip('\n') for line in open(str(n) + '_total.txt')]
    totals.append(total_file)

    # TOTAL CELL TEST
    total_file = [line.rstrip('\n') for line in open(str(n) + '_totalcell_test.txt')]
    total_cell_tests.append(total_file)

    # TYPE 1
    type1_file = [line.rstrip('\n') for line in open(str(n) + '_type1.txt')]
    type1s.append(type1_file)
    # TYPE 1R
    type1_r_file = [line.rstrip('\n') for line in open(str(n) + '_type1_r.txt')]
    type1_rs.append(type1_r_file)
    # TYPE 2
    type2_file = [line.rstrip('\n') for line in open(str(n) + '_type2.txt')]
    type2s.append(type2_file)
    # TYPE 2R
    type2_r_file = [line.rstrip('\n') for line in open(str(n) + '_type2_r.txt')]
    type2_rs.append(type2_r_file)
    # TYPE 3
    type3_file = [line.rstrip('\n') for line in open(str(n) + '_type3.txt')]
    type3s.append(type3_file)

    # OXYGEN
    oxygen_file = [line.rstrip('\n') for line in open(str(n) + '_oxygen_test.txt')]
    oxygens.append(oxygen_file)
    # CHEMO
    chemo_file = [line.rstrip('\n') for line in open(str(n) + '_chemo1.txt')]
    chemos.append(chemo_file)
    # CKINE
    ckine_file = [line.rstrip('\n') for line in open(str(n) + '_ckine.txt')]
    ckines.append(ckine_file)
    # CONTENTS
    contents_file = [line.rstrip('\n') for line in open(str(n) + '_data_test.txt')]
    contents.append(contents_file)

# Data gathered, now sum & write to files
activemac_file = open('output/activemac.txt', 'a')
caseation_file = open('output/caseation.txt', 'a')
chroninfectedmac_file = open('output/chroninfectedmac.txt', 'a')
infectedmac_file = open('output/infectedmac.txt', 'a')
intra_bac_file = open('output/intra_bac.txt', 'a')
restingmac_file = open('output/restingmac.txt', 'a')
total_file = open('output/total.txt', 'a')
total_cell_test_file = open('output/totalcell_test.txt', 'a')
type1_file = open('output/type1.txt', 'a')
type1_r_file = open('output/type1_r.txt', 'a')
type2_file = open('output/type2.txt', 'a')
type2_r_file = open('output/type2_r.txt', 'a')
type3_file = open('output/type3.txt', 'a')

for index in range(len(activemacs[0])):
    activemac_file.write(str(sum([int(activemacs[n][index]) for n in range(number_of_tiles)])))
    activemac_file.write('\n')
    caseation_file.write(str(sum([int(caseations[n][index]) for n in range(number_of_tiles)])))
    caseation_file.write('\n')
    chroninfectedmac_file.write(str(sum([int(chroninfectedmacs[n][index]) for n in range(number_of_tiles)])))
    chroninfectedmac_file.write('\n')
    infectedmac_file.write(str(sum([int(infectedmacs[n][index]) for n in range(number_of_tiles)])))
    infectedmac_file.write('\n')
    intra_bac_file.write(str(sum([int(intra_bacs[n][index]) for n in range(number_of_tiles)])))
    intra_bac_file.write('\n')
    restingmac_file.write(str(sum([int(restingmacs[n][index]) for n in range(number_of_tiles)])))
    restingmac_file.write('\n')
    total_file.write(str(sum([int(totals[n][index]) for n in range(number_of_tiles)])))
    total_file.write('\n')
    type1_file.write(str(sum([int(type1s[n][index]) for n in range(number_of_tiles)])))
    type1_file.write('\n')
    type1_r_file.write(str(sum([int(type1_rs[n][index]) for n in range(number_of_tiles)])))
    type1_r_file.write('\n')
    type2_file.write(str(sum([int(type2s[n][index]) for n in range(number_of_tiles)])))
    type2_file.write('\n')
    type2_r_file.write(str(sum([int(type2_rs[n][index]) for n in range(number_of_tiles)])))
    type2_r_file.write('\n')
    type3_file.write(str(sum([int(type3s[n][index]) for n in range(number_of_tiles)])))
    type3_file.write('\n')

for index in range(len(total_cell_tests[0])):
    total_cell_test_file.write(str(sum([int(total_cell_tests[n][index]) for n in range(number_of_tiles)])))
    total_cell_test_file.write('\n')

# Grid values
# TODO - assumes 2D 2x2 tile arrangement
x = 0
y = 0
counters = [0,] * number_of_tiles

oxygen_file = open('output/oxygen_test.txt', 'a')
chemo_file = open('output/chemo1.txt', 'a')
ckine_file = open('output/ckine.txt', 'a')
contents_file = open('output/data_test.txt', 'a')


for index in range(10000 * 200):

    if x < tile_shape[0] and y < tile_shape[1]:
        current_tile = 0
    elif x < tile_shape[0] and y >= tile_shape[1]:
        current_tile = 1
    elif x >= tile_shape[0] and y < tile_shape[1]:
        current_tile = 2
    elif x >= tile_shape[0] and y >= tile_shape[1]:
        current_tile = 3
    else:
        raise Exception('FAIL')

    oxygen_file.write(str(oxygens[current_tile][counters[current_tile]]))
    oxygen_file.write('\n')
    chemo_file.write(str(chemos[current_tile][counters[current_tile]]))
    chemo_file.write('\n')
    ckine_file.write(str(ckines[current_tile][counters[current_tile]]))
    ckine_file.write('\n')
    contents_file.write(str(contents[current_tile][counters[current_tile]]))
    contents_file.write('\n')

    counters[current_tile] += 1
    y += 1
    if y == 100:
        y = 0
        x += 1
        if x == 100:
            x = 0
