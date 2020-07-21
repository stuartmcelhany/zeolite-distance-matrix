import csv
import math

'''
ENTER INFORMATION HERE
'''
FILE = 'isomer-1.xyz'   # name of unit cell file
SIDE_LENGTH = 17.53     # side length of unit cell
ANGLE = 60              # angle of unit cell in degrees


'''
Creates primitive vectors for trigonal unit cell
'''
ANGLE = ANGLE * (math.pi/180)
PRIM_VEC_X = [SIDE_LENGTH, 0.0, 0.0] # x translation (x, y, z)
PRIM_VEC_Y = [SIDE_LENGTH*math.cos(ANGLE), SIDE_LENGTH*math.sin(ANGLE), 0.0] # y translation
cx = SIDE_LENGTH*math.cos(ANGLE)
cy = SIDE_LENGTH*(math.cos(ANGLE)-math.cos(ANGLE)**2)/math.sin(ANGLE)
PRIM_VEC_Z = [cx, cy, math.sqrt(SIDE_LENGTH**2-cx**2-cy**2)] # z translation


'''
Removes uncessary spaces from rows
'''
def clean_row_up(row):
    for value in list(row):
        if value == '':
            row.remove(value)
    return row


'''
Detects if there is a bond between Al or Si
'''
def get_bond(A, B):
    dist = math.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2+(B[2]-A[2])**2)
    if dist > 3.00 and dist < 3.60:
        return 1 # there is a bond
    else:
        return 0 # there is NOT a bond


'''
Creates node and edge graph in dictionary format
'''
def generate_graph(molecule):
    graph = {}
    for row in molecule:
        graph[row[5]] = {} # creates key
        coordsA = [float(row[1]), float(row[2]), float(row[3])]
        for row_temp in molecule:
            coordsB = [float(row_temp[1]), float(row_temp[2]), float(row_temp[3])]
            if get_bond(coordsA, coordsB):
                graph[row[5]][row_temp[5]] = 1
    return graph


'''
Run the Dijkstra algorithm to find shortest path
Note that each weight is simply 1, indicating there exists a bond between two atoms
'''
def dijkstra(graph,start,goal):
    shortest_distance = {}
    predecessor = {}
    unseenNodes = graph.copy()
    infinity = 9999999
    path = []
    for node in unseenNodes:
        shortest_distance[node] = infinity
    shortest_distance[start] = 0
 
    while unseenNodes:
        minNode = None
        for node in unseenNodes:
            if minNode is None:
                minNode = node
            elif shortest_distance[node] < shortest_distance[minNode]:
                minNode = node
 
        for childNode, weight in graph[minNode].items():
            if weight + shortest_distance[minNode] < shortest_distance[childNode]:
                shortest_distance[childNode] = weight + shortest_distance[minNode]
                predecessor[childNode] = minNode
        unseenNodes.pop(minNode)
 
    currentNode = goal
    while currentNode != start:
        try:
            path.insert(0,currentNode)
            currentNode = predecessor[currentNode]
        except KeyError:
            return infinity
    path.insert(0,start)
    if shortest_distance[goal] != infinity:
        #print('Shortest distance is ' + str(shortest_distance[goal]))
        #print('And the path is ' + str(path))
        return shortest_distance[goal] - 1


'''
3x3x3 Cloner that removes all other atoms except Al and Si
This builds the molecule file in periodic conditions
'''
mol = []
with open(FILE) as File:
    reader = csv.reader(File, delimiter = ' ', lineterminator = '\n')
    next(reader) # skips atom count line
    next(reader) # skips comment line per .xyz convention
    original_index = 1
    new_index = 0
    for row in reader:
        row = clean_row_up(row)
        if row[0] == 'Al' or row[0] == 'Si':
            new_index += 1
            row_entry = [""]*6
            row_entry[0] = row[0] #get atom
            row_entry[1] = "{:.4f}".format(float(row[1])) # get x
            row_entry[2] = "{:.4f}".format(float(row[2])) # get y
            row_entry[3] = "{:.4f}".format(float(row[3])) # get z
            row_entry[4] = original_index
            row_entry[5] = new_index
            mol.append(row_entry)
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        if i==1 and j==1 and k==1:
                            continue
                        else:
                            new_index += 1
                            row_entry = [""]*6
                            row_entry[0] = row[0] #get atom
                            row_entry[1] = "{:.4f}".format(float(row[1]) + float((i-1)*PRIM_VEC_X[0]) + float((j-1)*PRIM_VEC_Y[0]) + float((k-1)*PRIM_VEC_Z[0])) # get x
                            row_entry[2] = "{:.4f}".format(float(row[2]) + float((i-1)*PRIM_VEC_X[1]) + float((j-1)*PRIM_VEC_Y[1]) + float((k-1)*PRIM_VEC_Z[1])) # get y
                            row_entry[3] = "{:.4f}".format(float(row[3]) + float((i-1)*PRIM_VEC_X[2]) + float((j-1)*PRIM_VEC_Y[2]) + float((k-1)*PRIM_VEC_Z[2])) # get z
                            row_entry[4] = original_index
                            row_entry[5] = new_index
                            mol.append(row_entry)
        original_index += 1


'''
Generates inter-tetrahedron distance matrix
'''
increaseCount = 0 # used for loading bar
Al_count = 0 # used for loading bar
for i, row in enumerate(mol): # used for loading bar
    if row[0] == 'Al' and (i % 27) == 0: # used for loading bar
        Al_count += 1 # used for loading bar
count = 100/Al_count # used for loading bar
graph = generate_graph(mol)
dist_array = [['-']]
line_indexA = 0
Al_index = 0
for row in mol:
    line_indexB = 0
    if row[0] == 'Al' and (line_indexA %  27) == 0:
        Al_index += 1
        print('Loading: '+str(int(increaseCount))+'%',end='\r')
        increaseCount += count
        dist_array[0].append('Al' + str(row[4]))
        dist_array.append(['Al' + str(row[4])])
        for newRow in mol:
            temp_dist = []
            if line_indexA == line_indexB:
                dist_array[Al_index].append('-')
            elif newRow[0] == 'Al' and (line_indexB % 27) == 0:
                for i in range(27):
                    temp_dist.append(dijkstra(graph, line_indexA+1, line_indexB+i+1)) # indices have 1 added to them to match node naming format
                dist_array[Al_index].append(min(temp_dist))
            line_indexB += 1
    line_indexA += 1


'''
Write to output csv file
'''
FILE = FILE.replace('.xyz', '')
with open(FILE + '-matrix.csv', 'w') as File:
    writer = csv.writer(File, delimiter = ',', lineterminator = '\n') # ',' delimiter is used for columns in excel
    writer.writerows(dist_array)