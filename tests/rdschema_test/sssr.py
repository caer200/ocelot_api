from collections import deque

"""
got this from https://gist.github.com/zachcp/c39afdc2353f2412ae50

unfortunately this does not work for rings connected by a path, see the case of four_six_plus_three
"""

pentagon = [{"key":"A", "connections":["B","E"]},
   {"key":"B", "connections": ["A", "C"]},
   {"key":"C", "connections": ["B", "D"]},
   {"key":"D", "connections": ["C", "E"]},
   {"key":"E", "connections": ["D", "A"]}]

hexagon = [{"key":"A","connections":["B", "E"]},
   {"key":"B","connections":["A", "C"]},
   {"key":"C","connections":["B", "D"]},
   {"key":"D","connections":["C", "E"]},
   {"key":"E","connections":["D", "A"]}]

four_six = [
{"key":"A","connections":["B", "H"]},
{"key":"B","connections":["A", "C"]},
{"key":"C","connections":["B", "D", "H"]},
{"key":"D","connections":["C", "E"]},
{"key":"E","connections":["D", "F"]},
{"key":"F","connections":["E", "G"]},
{"key":"G","connections":["F", "H"]},
{"key":"H","connections":["A", "C", "G"]}
]

three_five = [
{"key":"A","connections":["B", "F"]},
{"key":"B","connections":["A", "C", "F"]},
{"key":"C","connections":["B", "D"]},
{"key":"D","connections":["C", "E"]},
{"key":"E","connections":["D", "F"]},
{"key":"F","connections":["E", "A", "B"]}
]

four_six_plus = [
{"key":"A","connections":["B", "H", "I"]},
{"key":"B","connections":["A", "C"]},
{"key":"C","connections":["B", "D", "H"]},
{"key":"D","connections":["C", "E"]},
{"key":"E","connections":["D", "F"]},
{"key":"F","connections":["E", "G"]},
{"key":"G","connections":["F", "H"]},
{"key":"H","connections":["A", "C", "G"]},
{"key":"I","connections":["A","J"]},
{"key":"J","connections":["I"]}
]

four_six_plus_three = [
    {"key":"A","connections":["B", "H", "I"]},
    {"key":"B","connections":["A", "C"]},
    {"key":"C","connections":["B", "D", "H"]},
    {"key":"D","connections":["C", "E"]},
    {"key":"E","connections":["D", "F"]},
    {"key":"F","connections":["E", "G"]},
    {"key":"G","connections":["F", "H"]},
    {"key":"H","connections":["A", "C", "G"]},
    {"key":"I","connections":["A","J"]},
    {"key":"J","connections":["I", "K", "M"]},
    {"key":"K","connections":["J", "L"]},
    {"key": "L", "connections": ["K", "M"]},
    {"key": "M", "connections": ["L", "J"]}

]

def getring(Molecule, rootnode):
    # connections is just the table of connections in the molecule
    connections = { mol['key']: mol['connections'] for mol in Molecule}
    # ring members must be a member of the ring
    assert len(connections[rootnode])>1
    # paths is a dictionary of paths starting form the start node. It is initialized
    # to zero and will be used to updated as each member of the queue is analyzed
    paths = {mol['key']: set() for mol in Molecule}
    
    # atomqueue will keep track of atoms that have been visited and need to be processed
    atomqueue = deque()
    
    # visited is the atom keys that have been seen already
    visited = set()
    
    #when a ring is found, return the ring
    ring = None
    
    # initialize the root node
    atomqueue.append(rootnode)
    paths[rootnode].add(rootnode)
    visited.add(rootnode)

    while not ring:
        print ("atomqueue: {}".format(atomqueue))
        print ("visited: {}".format(visited))
        print ("paths: {}".format(paths))
        
        try:
            atomcode = atomqueue.popleft()
            #print "atomcode: {}".format(atomcode)

            connectedatoms = [a for a in connections[atomcode] if a not in visited]
            #print "atoms connected: {}".format(connectedatoms)
            for atom in connectedatoms:
                print ("atom: {}".format(atom))
                #check for rings and return if a ring is found
                currentpath = paths[atomcode]
                nextpath  = paths[atom]
                intersect = currentpath.intersection(nextpath)
                #print "current path: {}".format(currentpath)
                #print "next path: {}".format(currentpath)
                if len(nextpath) > 0 and len(intersect) == 1:
                    return currentpath.union(nextpath)
    
                # update path only if the target path is empty 
                # to avoid incorporating longer paths
                if len(nextpath) == 0:
                    paths[atom] = currentpath.copy()
                    paths[atom].add(atom)
            
                #add atoms to queue
                atomqueue.append(atom)
                
            #update visited
            visited.add(atomcode)
        except:
            print("Problem Here. You should never get here!")

def getrings(Molecule):
    """
    find atoms that are connected by only a single bond,
    eliminate them and the connections to them by other atoms until
    there are only atoms that have at least two connections
    then find the rings on them
    """
    
    # copy the Molecule for local editing
    tempmol = list(Molecule)
    print (tempmol)
    # find the atoms that are initially less than 2 bond connections
    singles_index = [i for i,a in enumerate(tempmol) if len(a['connections']) < 2]
    
    
    # eliminate all atoms that are not part of a ring
    while len(singles_index) > 0:
        for idx in singles_index:
            key = tempmol[idx]['key']
            connections = tempmol[idx]['connections']
            
            # remove connections 
            for con in connections:
                #print "index: {}\nkey:{}\nconnection:{}".format(idx,key,con)
                index = [i for i,x in enumerate(tempmol) if x['key'] == con][0]
                tempmol[index]['connections'].remove(key)
 
            # remove atom
            tempmol.pop(idx) 
            print (tempmol)
            
            # reset singles index
            singles_index = [i for i,a in enumerate(tempmol) if len(a['connections']) < 2]
            print (len(singles_index))
        
    # calculate ring on the trimmed molceule
    rings = set()
    for atom in tempmol:
        if len(atom['connections']) > 1:
            ring = frozenset(getring(Molecule, atom["key"]))
            print (atom, ring)
            if ring not in rings:
                rings.add(ring)
    return rings


# getring(pentagon, 'B') #expect six
rings = getrings(four_six_plus_three)
print(rings)
# getring(hexagon, 'B')  # expect five
# getring(four_six, 'A') #expect four yes
# getring(four_six, 'B') #expect four yes
# getring(four_six, 'C') #expect four yes
# getring(four_six, 'D') #expect six yes
# getring(four_six, 'E') #expect six yes
# getring(four_six, 'F') #expect six yes
# getring(four_six, 'G') #expect siz yes
# getring(four_six, 'H') #expect four yes
# getring(four_six_plus, 'I') #expect assertion for length
# getrings(four_six) # correctly gets rings
# getrings(three_five) # incorrect for C,E,F
# getrings(four_six_plus) # should trim atoms and correctly get rings
