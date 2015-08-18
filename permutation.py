class Permutation(object):
    def __init__(self, sequence):
        self.permutation = sequence

    def getPermuted(self, sequence):
        return [sequence[i] for i in self.permutation]

    def getInversePermuted(self, sequence):
        return [sequence[self.permutation.index(i)] for i in range(len(self.permutation))]
