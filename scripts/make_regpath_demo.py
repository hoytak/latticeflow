#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV, calculate2dTVPath
from matplotlib.pylab import imread, figure, show
from matplotlib import collections
from itertools import product
from time import clock

lambda_bounds = [0.1, .25]
lineweight = 0.05


plot_type = "levels2"

#image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/branches.png"
#image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/sanity.png"
#image_file = "benchmarks/images/branches-small.png"
#image_file = "benchmarks/images/mona_lisa.png"
image_file = "benchmarks/images/Bodiam-castle.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())[::2, ::2]

X -= X.mean()
X /= X.std()

# X = np.linspace(0,1.25,6).reshape( (3,2) )

# Xtv0 = calculate2dTV(X, 0)

# assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

tv_path_lambdas = np.linspace(lambda_bounds[0], lambda_bounds[1], 100)

time_start = clock()
Xtv_path = calculate2dTVPath(X, tv_path_lambdas)
tv_path_time = clock() - time_start

print "Done calculating Regpath Version; time = %f" % tv_path_time

tv_lambdas = np.linspace(lambda_bounds[0], lambda_bounds[1], 100)

def plotRegPath(a, lambdas, Xtv, transform):

    xc, yc = 0, 0

    col = []
    
    ymin, ymax = 0,0



    for xi, yi in product(range(0, Xtv.shape[1]), range(0, Xtv.shape[2])):


        y = Xtv[:,xi, yi]
        
        if transform:
            y /= lambdas

        y += np.random.normal(0,0.001)
        
        col.append( np.vstack([lambdas, y]).T )
        
        ymin = min(ymin, y.min())
        ymax = max(ymax, y.max())

    a.add_collection(collections.LineCollection(col, antialiased = True, linewidths=lineweight))
    a.legend()
    a.set_xlim([lambdas.min(),lambdas.max()])
    a.set_ylim([ymin, ymax])

f = figure()
a = f.add_subplot(111)
plotRegPath(a, tv_path_lambdas, Xtv_path, False)
a.set_title(r"Raw Reg Path on $\lambda \in \Icc{%0.2f, %0.2f}$" % (tv_path_lambdas.min(), tv_path_lambdas.max()))

f = figure()
a = f.add_subplot(111)
plotRegPath(a, tv_path_lambdas, Xtv_path, True)
a.set_title(r"Reg Path on $\lambda \in \Icc{%0.2f, %0.2f}$" % (tv_path_lambdas.min(), tv_path_lambdas.max()))

show()


# Get the first round
Xtv_1 = np.empty( (tv_lambdas.size, X.shape[0], X.shape[1]) )

time_start = clock()
for i, lm in enumerate(tv_lambdas):
    print "Calculating %d/%d (lambda = %1.5f)" % ((i + 1), len(tv_lambdas), lm)
    Xtv_1[i,:,:] = calculate2dTV(X, lm)

tv_round_time = clock() - time_start

print "Done calculating indivdual models; time took %f." % tv_round_time
