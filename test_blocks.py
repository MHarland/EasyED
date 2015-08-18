from numpy import array

from blocks import BlockMatrix, BlockState, bdot

x = BlockMatrix([2,2,1])
a = array([[1,2],[3,4]])
b = array([[1,5],[5,4]])
c = array([[10]])
abc = [a, b, c]
for i, block in enumerate(x):
    x.datablocks[i] = abc[i]

x.show()
print

y = bdot(x, x)
y.show()
print

z = BlockState([2,2,1])
for i, block in enumerate(z):
    z.datablocks[i] = abc[i][0, :]
z.show()
print

y = bdot(x,z)
y.show()
print

x.setFullMatrixEntry(2,3,7)
x.setFullMatrixEntry(1,1,7)
x.setFullMatrixEntry(0,1,7)
x.setFullMatrixEntry(4,4,7)
x.setFullMatrixEntry(3,3,7)
x.setFullMatrixEntry(3,2,7)
x.setFullMatrixEntry(0,0,7)
x.setFullMatrixEntry(1,0,7)
x.setFullMatrixEntry(2,2,7)
x.show()

x.setFullMatrixEntry(0,2,7) #error
