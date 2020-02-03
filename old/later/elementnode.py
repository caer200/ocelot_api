from pymatgen.core.periodic_table import Element

class ElementNode:
    """
    this is just a tuple of (name, index), it should be immutable
    """

    __slots__ = ('name', 'index')

    def __init__(self, name, index):
        if not isinstance(name, str):
            raise ValueError("'name' must be a string")
        if not isinstance(index, int):
            raise ValueError("'index' must be an int")
        super(ElementNode, self).__setattr__('name', name)
        super(ElementNode, self).__setattr__('index', int(index))

    def __setattr__(self, name, value):
        raise AttributeError('Persons cannot be modified')

    def __repr__(self):
        return "{}: {}, {}".format(self.__class__.__name__, self.name, self.index)

    def __hash__(self):
        return hash((self.name, self.index))

    def __eq__(self, other):
        return (self.name, self.index) == (other.name, other.index)

    @property
    def element(self):
        try:
            return Element(self.name)
        except ValueError:
            raise ValueError("{} is not a valid element name".format(self.name))




