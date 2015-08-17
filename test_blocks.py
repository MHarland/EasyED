from numpy import array

from blocks import BlockMatrix, BlockState, bdot

x = BlockMatrix([2,2,1])
a = array([[1,2],[3,4]])
b = array([[1,5],[5,4]])
c = array([[10]])
abc = [a, b, c]
for i, block in enumerate(x):
    x.datablocks[i] = abc[i]

for block in x:
    print block
print

y = bdot(x, x)
for block in y:
    print block
print

z = BlockState([2,2,1])
for i, block in enumerate(z):
    z.datablocks[i] = abc[i][0, :]
for block in z:
    print block
print

y = bdot(x,z)
for block in y:
    print block
print

x.setFullMatrixEntry(2,3,7)
for block in x:
    print block
