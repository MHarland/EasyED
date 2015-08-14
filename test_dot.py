from util import dot
from numpy import array

a = array([1,2])
b = array([[1,2,],[3,4]])
c = array([1,2])

print dot(a,dot(b,c))
print dot(a,b,c)
print dot(dot(a,b),c)

print dot(b,c)
