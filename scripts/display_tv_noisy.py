#!/use/bin/env python

import numpy as np
import sys
from pylatticeflow import calculate2dTV
from matplotlib.pylab import imread, figure, show

lm = float(sys.argv[1])
#image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/trees-cut.png"
#image_file = "benchmarks/images/marmot.png"
#image_file = "benchmarks/images/mona_lisa.png"
#image_file = "benchmarks/images/branches.png"
#image_file = "benchmarks/images/Bodiam-castle.png"
image_file = "benchmarks/images/ct-brain.png"


Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

print "Image file %s loaded." % image_file

X = (Xo.mean(axis=2) / Xo.max())

X -= X.mean()
X /= X.std()

X_max = X.max()
X_min = X.min()

X += np.random.normal(0, 0.25, size = X.shape)

Xtv = calculate2dTV(X, lm)

Xtv -= X.mean()
Xtv /= X.std()

# Xtv[Xtv > X_max] = X_max
# Xtv[Xtv < X_min] = X_min

# X += X.min()
# X /= X.max()

# Xtv += X.min()
# Xtv /= X.max()


f = figure()
a = f.add_subplot(111)
a.set_title("Original Image.")
a.imshow(X, interpolation='nearest', cmap="gray")

f = figure()
a = f.add_subplot(111)
a.set_title("Total Variation Solution, $\lambda=%0.2f$." % lm)
a.imshow(Xtv,  interpolation='nearest', cmap="gray")

show()

