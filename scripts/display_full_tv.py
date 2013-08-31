#!/use/bin/env python

import numpy as np
import sys
from pylatticeflow import calculate2dTVImage
from matplotlib.pylab import imread, figure, show

lm = float(sys.argv[1])
#image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/trees.png"
#image_file = "benchmarks/images/marmot.png"
image_file = "benchmarks/images/mona_lisa.png"
#image_file = "benchmarks/images/branches.png"

X = imread(image_file)

if not X.size:
    raise IOError("Error loading image %s." % image_file)

print "Image file %s loaded." % image_file

# X -= X.mean()
# X /= X.std()

Xtv = calculate2dTVImage(X, lm)

# Xtv -= X.mean()
# Xtv /= X.std()

f = figure()
a = f.add_subplot(111)
a.set_title("Original Image.")
a.imshow(X, interpolation='nearest')

f = figure()
a = f.add_subplot(111)
a.set_title("Total Variation Levels.")
a.imshow(Xtv,  interpolation='nearest')


show()

