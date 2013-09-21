#!/use/bin/env python

import numpy as np
import sys
from itertools import product
from pylatticeflow import partiallyObservedRegression
from matplotlib.pylab import imread, figure, show

# lm = float(sys.argv[1])
#image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/trees-cut.png"
#image_file = "benchmarks/images/marmot.png"
#image_file = "benchmarks/images/mona_lisa.png"
#image_file = "benchmarks/images/branches.png"
#image_file = "benchmarks/images/Bodiam-castle.png"

reg_p = float(sys.argv[1])
n_iterations = float(sys.argv[2])

n_slices = 8
std_thresh = 0

image_file = "benchmarks/images/ct-brain.png"
np.random.seed(0)

max_n_opservations = 1000

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

print "Image file %s loaded." % image_file

X = (Xo.mean(axis=2) / Xo.max())[::5, ::5]

X -= X.mean()
X /= X.std()


y_idx = set()

def add_set(xl, xh, yl, yh):
    for i, j in product(xrange(xl, xh), xrange(yl,yh)):
        if 0 <= i < X.shape[0] and 0 <= j < X.shape[1]:
            y_idx.add( (i, j) )

#b = int(min(X.shape[0], X.shape[1]) )

# add_set(0, b, 0, X.shape[1])
# add_set(X.shape[0]-b, X.shape[0], 0, X.shape[1])
# add_set(0, X.shape[0], 0, b)
# add_set(0, X.shape[0], X.shape[1] - b, X.shape[1])

def add_block(x, y):
    if not (-1 <= x <= 1 and -1 <= y <= 1):
        return False
    
    ix = int((x + 1.0) *0.5 * X.shape[0])
    iy = int((y + 1.0) *0.5 * X.shape[1])

    samples = []

    for i, j in product(xrange(ix - 2, ix + 3), xrange(iy - 2, iy + 3)):
        if 0 <= i < X.shape[0] and 0 <= j < X.shape[1]:
            samples.append( X[i,j] )

    assert samples

    if np.array(samples).std() > std_thresh:
        add_set(ix - 1, ix + 2, iy - 1, iy + 2)

    return True

for i in xrange(n_slices):

    angle = i * 3.14159 / n_slices
    xv = np.cos(angle)
    yv = np.sin(angle)

    for t in np.linspace(0,2, 1000):
        add_block(t*xv, t*yv)
        add_block(-t*xv, -t*yv)
    
y_idx = sorted(list(y_idx))

n_observations = len(y_idx)

if n_observations > max_n_opservations:
    n_observations = max_n_opservations
    elements = np.random.choice(range(len(y_idx)), n_observations, replace = False)
    y_idx = [y_idx[i] for i in elements]

value_matrix = np.zeros(X.shape)

for i, j in y_idx:
    value_matrix[i,j] = X[i,j]

f = figure()
a = f.add_subplot(111)
a.set_title("Original Image.")
a.imshow(X, interpolation='nearest', cmap="gray")

f = figure()
a = f.add_subplot(111)
a.set_title("Inpainting Solution")
a.imshow(value_matrix,  interpolation='nearest', cmap="gray")

print "Starting regularization." 

Xtr, path = partiallyObservedRegression(value_matrix, reg_p, n_iterations, 1)

f = figure()
a = f.add_subplot(111)
a.set_title("FISTA solution") 
a.imshow(Xtr,  interpolation='nearest', cmap="gray")

show()


