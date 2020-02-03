def check_continuity(my_list, delta=1):
    return all(a + delta == b for a, b in zip(my_list, my_list[1:]))


class NodeCollection:

    def __init__(self, nodes):
        """
        only nodes with unique ids will be added into the collection
        """
        self.nodes = [n for n in self.unique_by_index(nodes)]
        self.nodes.sort(key=lambda x: x.index)
        self.nodes = tuple(self.nodes)
        self.indices = [n.index for n in self.nodes]

    @property
    def idcontinuous(self):
        """
        check whether the indices of nodes are continuous
        """
        return check_continuity(self.nodes)

    def __contains__(self, node):
        for n in self.nodes:
            if n == node:
                return True
        return False

    @staticmethod
    def unique_by_index(nodes):
        seen = set()
        for node in nodes:
            if not node.index in seen:
                seen.add(node.index)
                yield node

    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        return self.nodes.__iter__()

    def __getitem__(self, index):
        for n in self.nodes:
            if n.index == index:
                return n
        raise IndexError('index {} is not found in this collection'.format(index))

    def __repr__(self):
        outs = [self.__class__.__name__ + ': ']
        for s in self.nodes:
            outs.append(s.__repr__())
        return '\n'.join(outs)
