from numpy import *
from pylab import *

x = loadtxt('fort.45')
y = x[:,3]
m = size(y)
n = int(sqrt(m))
y = x[:,3].reshape(n,n)
figure(1)
imshow(y,extent=[-3,3,-3,3])
colorbar()
show()

figure(2)
x = loadtxt('fort.46')
y = x[:,3]
m = size(y)
n = int(sqrt(m))
y = x[:,3].reshape(n,n)
imshow(y,extent=[-3,3,-3,3])
colorbar()
show()
