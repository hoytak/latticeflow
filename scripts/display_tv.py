#!/use/bin/env python

import numpy as np
import sys
from pylatticeflow import calculate2dTV
from matplotlib.pylab import imread, figure, show

lm = float(sys.argv[1])
image_file = "benchmarks/images/branches.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())

Xtv = calculate2dTV(X, lm)

f = figure()
a = f.add_subplot(121)
a.imshow(X, vmin=0, vmax=1)

a = f.add_subplot(122)
a.imshow(Xtv, vmin=0, vmax=1)

show()

