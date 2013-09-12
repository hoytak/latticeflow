#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV, calculate2dTVPath
from matplotlib.pylab import imread, figure, show
from matplotlib import collections
from itertools import product

plot_type = "levels2"

#image_file = "benchmarks/images/truffles.png"
image_file = "benchmarks/images/truffles-small.png"
# image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/branches-small.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())

X -= X.mean()
X /= X.std()

# X = np.linspace(0,1.25,6).reshape( (3,2) )

# Xtv0 = calculate2dTV(X, 0)

# assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

lambdas = np.linspace(0.001, 0.1, 500)

Xtv_2 = calculate2dTVPath(X, lambdas)
print "Done calculating Regpath Version."

# # Get the first round
# Xtv_1 = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )

# for i, lm in enumerate(lambdas):
#     print "Calculating %d/%d (lambda = %1.5f)" % ((i + 1), len(lambdas), lm)
#     Xtv_1[i,:,:] = calculate2dTV(X, lm)

# print "Done calculating indivdual models."


def plotRegPath(a, Xtv):

    xc, yc = 0, 0

    col = []
    
    ymin, ymax = 0,0

    for xi, yi in product(range(0, Xtv.shape[1]), range(0, Xtv.shape[2])):
        
        y = Xtv[:,xi, yi] / lambdas

        col.append( np.vstack([lambdas, y]).T )
        
        ymin = min(ymin, y.min())
        ymax = max(ymax, y.max())

    a.add_collection(collections.LineCollection(col, antialiased = True, linewidths=0.5))
    a.legend()
    a.set_xlim([lambdas.min(),lambdas.max()])
    a.set_ylim([ymin, ymax])

f = figure()

plotRegPath(f.add_subplot(111), Xtv_2 / Xtv_2.std())

show()
