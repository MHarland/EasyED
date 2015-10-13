from permutation import Permutation

p = Permutation([3,0,1,2])
x = p.getPermuted(range(4))
print x
range4 = p.getInversePermuted(x)
print range4
