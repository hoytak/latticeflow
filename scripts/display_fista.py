#!/use/bin/env python

import numpy as np
import sys
from itertools import product
from pylatticeflow import regularizedRegression
from matplotlib.pylab import imread, figure, show

# lm = float(sys.argv[1])
#image_file = "benchmarks/images/sanity.png"
image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/trees-cut.png"
#image_file = "benchmarks/images/marmot.png"
#image_file = "benchmarks/images/mona_lisa.png"
#image_file = "benchmarks/images/branches.png"
#image_file = "benchmarks/images/Bodiam-castle.png"

np.random.seed(0)

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

print "Image file %s loaded." % image_file

X = (Xo.mean(axis=2) / Xo.max())

X -= X.mean()
X /= X.std()

n_observations = 100

y_idx = set()

def add_set(xl, xh, yl, yh):
    for i, j in product(xrange(xl, xh), xrange(yl,yh)):
        y_idx.add( (i, j) )


b = int(min(X.shape[0], X.shape[1]) / 4)

add_set(0, b, 0, X.shape[1])
add_set(X.shape[0]-b, X.shape[0], 0, X.shape[1])
add_set(0, X.shape[0], 0, b)
add_set(0, X.shape[0], X.shape[1] - b, X.shape[1])

y_idx = sorted(list(y_idx))

n_observations = len(y_idx)

A = np.zeros( X.shape + (n_observations,) ) 
y = np.empty(n_observations)
Obs = np.zeros(X.shape)

for i in xrange(n_observations):
    xi, yi = y_idx[i]
    A[xi,yi,i] = 1
    y[i] = Obs[xi,yi] = X[xi,yi]

print "Starting regularization." 

Xtr, path = regularizedRegression(A, y, 1, 100, 1)

f = figure()
a = f.add_subplot(111)
a.set_title("Original Image.")
a.imshow(X, interpolation='nearest', cmap="gray")

f = figure()
a = f.add_subplot(111)
a.set_title("Inpainting Solution")
a.imshow(Obs,  interpolation='nearest', cmap="gray")

f = figure()
a = f.add_subplot(111)
a.set_title("FISTA solution") 
a.imshow(Xtr,  interpolation='nearest', cmap="gray")

show()

