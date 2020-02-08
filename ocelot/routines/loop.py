"""
as we are using networkx this can be deprecated...
"""

class Loopsearcher:

    def __init__(self, nbmap):
        """
        giving a connection table, find all possible loops with a certain loop size

        nbmap[i] does not contain i itself

        useful to consider smallest set of smallest rings (sssr) problem

        check 10.1073pnas.0813040106 for a better solution

        :param nbmap: connection table
        """
        self.nbmap = nbmap
        self.edges = self.generate_edges()
        self.nodes = list(self.nbmap.keys())

    def expand(self, path_set):
        """
        the path will never intersect itself
        """
        new_path_set = []

        if len(path_set) == 1 and len(path_set[0]) == 1:
            start = path_set[0][0]
            for n in self.nbmap[start]:
                if n != start:
                    new_path_set.append([start, n])
            return new_path_set

        for path in path_set:
            for n in self.nbmap[path[-1]]:
                if n not in path:
                    new_path_set.append(path + [n])

        if len(new_path_set) == 0:
            return path_set
        else:
            return new_path_set

    def alex_method(self, loopsize):
        """
        I figured this out but I'm sure I'm not the first

        :param int loopsize: ring size to look for
        :return: a list of index
        """
        all_idx = list(self.nbmap.keys())
        loop_found = []
        visited = []
        while len(visited) != len(all_idx):
            unvisited = [i for i in all_idx if i not in visited]
            start = unvisited[0]
            path_set = [[start]]
            for i in range(loopsize - 1):
                path_set = self.expand(path_set)
            for path in path_set:
                if path[0] in self.nbmap[path[-1]]:
                    if set(path) not in [set(loop) for loop in loop_found]:
                        loop_found.append(path)
                        visited += [p for p in path if p not in visited]
                    break
            if start not in visited:
                visited += [start]

        true_loop_found = []
        for idxlst in loop_found:
            if len(idxlst) == loopsize:
                true_loop_found.append(idxlst)
        return true_loop_found

    def sssr_alex(self, size_min, size_max):
        loops = []
        for size in range(size_min, size_max+1):
            loops += self.alex_method(size)
        basis_loops = []
        basis_loops_set = []
        for loop1 in loops:
            loop1fused = False
            for loop2 in loops:
                if set(loop1).issuperset(set(loop2)) and len(loop1) > len(loop2):
                    loop1fused = True
                    break
            if not loop1fused and set(loop1) not in basis_loops_set:
                    basis_loops.append(loop1)
                    basis_loops_set.append(set(loop1))

        # # there should be an easier way for this...
        # # https://stackoverflow.com/questions/32071425/
        # # I cannot find a gaussian elimination implementation based on xor
        # edgematrix = np.zeros((len(loops), len(self.edges)))
        # for i in range(len(loops)):
        #     l_edges = self.loop2edges(loops[i])
        #     for j in range(len(self.edges)):
        #         if self.edges[j] in l_edges:
        #             edgematrix[i, j] = 1
        # pl, u = lu(edgematrix, permute_l=True)
        # basis_loops = [loops[i] for i in np.where(u.any(axis=1))[0]]
        return basis_loops


    @staticmethod
    def loop2edges(loop):
        """
        notice here loop is sth like [1, 2, 3], where it is implied that 1 is connected to 3

        :param loop:
        :return: e.g. [{1, 2}, {2, 3}, {1, 3}]
        """
        path = loop + [loop[0]]
        edges = [set(pair) for pair in zip(path[:-1], path[1:])]
        return edges

    def generate_edges(self):
        graph = self.nbmap
        edges = []
        for node in graph:
            for neighbour in graph[node]:
                edges.append({node, neighbour})
        return edges
