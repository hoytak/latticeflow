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
image_file = "benchmarks/images/Bodiam-castle.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

print "Image file %s loaded." % image_file

X = (Xo.mean(axis=2) / Xo.max())

X -= X.mean()
X /= X.std()

Xtv = calculate2dTV(X, lm)

Xtv -= X.mean()
Xtv /= X.std()

# X += X.min()
# X /= X.max()

# Xtv += X.min()
# Xtv /= X.max()

Xnew = -( (Xtv[1:,:-1] - Xtv[:-1,:-1])**2 + (Xtv[:-1,1:] - Xtv[:-1,:-1])**2)**(0.25)



f = figure()
a = f.add_subplot(111)
a.set_title("Total Variation Edges, $\lambda=%0.2f$." % lm)
a.imshow(Xnew,  interpolation='spline16', cmap="gray")

show()
