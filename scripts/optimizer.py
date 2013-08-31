#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTVPath
from matplotlib.pylab import imread

plot_type = "levels2"

image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/truffles-small.png"
# image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/branches-small.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())[::2, ::2]

X -= X.mean()
X /= X.std()

# X = np.linspace(0,1.25,6).reshape( (3,2) )

# Xtv0 = calculate2dTV(X, 0)

# assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

lambdas = np.linspace(0.5, 0.78, 50)

calculate2dTVPath(X, lambdas)

