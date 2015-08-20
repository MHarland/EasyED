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

x[2,3]=7
x[1,1]=7
x[0,1]=7
x[4,4]=7
x[3,3]=7
x[3,2]=7
x[0,0]=7
x[1,0]=7
x[2,2]=7
x.show()
print x[2,2]
print
print x.getBlock(1)
print
print 'state:'
y = BlockState([2,2,1])
y[0] = 1
y[1] = 1
y[2] = 1
y[3] = 5
y[4] = 1
y.show()
print y[3]
print
print 'z from bdot:'
x.show()
print
y.show()
print
z = x.bdot(y)
z.show()
print
z = bdot(x,y)
z.show()
print 'a vector as result of leftmul:'
z = bdot(y,x)
z.show()
print 'for matrix:'
z = bdot(x,x)
z.show()
print 'error on purpose:'
x[0,2]=7 #error
