class Loopsearcher:

    def __init__(self, nbmap):
        """
        smallest set of smallest rings
        giving a connection table, find all possible loops with a certain loop size
        nbmap[i] does not contain i itself
        check 10.1073pnas.0813040106 for a better solution
        :param nbmap: connection table
        """
        self.nbmap = nbmap

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
        :param loopsize:
        :return:
        """
        all_idx = list(range(len(self.nbmap)))
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
        return loop_found

