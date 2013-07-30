#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV, calculate2dTVPath
from matplotlib.pylab import imread, figure, show
from matplotlib import collections
from itertools import product

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

lambdas = np.linspace(0.1, 1, 50)

Xtv_2 = calculate2dTVPath(X, lambdas)
print "Done calculating Regpath Version."

# Get the first round
Xtv_1 = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )

for i, lm in enumerate(lambdas):
    print "Calculating %d/%d (lambda = %1.5f)" % ((i + 1), len(lambdas), lm)
    Xtv_1[i,:,:] = calculate2dTV(X, lm)

print "Done calculating indivdual models."


def plotRegPath(a, Xtv):

    xc, yc = 0, 0

    col = []
    
    ymin, ymax = 0,0

    for xi, yi in product(range(0, Xtv.shape[1]), range(0, Xtv.shape[2])):
        
        y = Xtv[:,xi, yi]

        col.append( np.vstack([lambdas, y]).T )
        
        ymin = min(ymin, y.min())
        ymax = max(ymax, y.max())

    a.add_collection(collections.LineCollection(col, antialiased = True, linewidths=0.5))
    a.legend()
    a.set_xlim([lambdas.min(),lambdas.max()])
    a.set_ylim([ymin, ymax])

f = figure()

plotRegPath(f.add_subplot(411), Xtv_1)

plotRegPath(f.add_subplot(412), Xtv_2)

for i in range(Xtv_1.shape[0]):
    Xtv_1[i,:,:] -= Xtv_1[i,:,:].mean()
Xtv_1 /= Xtv_1.std()    

for i in range(Xtv_2.shape[0]):
    Xtv_2[i,:,:] -= Xtv_2[i,:,:].mean()
Xtv_2 /= Xtv_2.std()    

print "Xtv_1.std() = ", Xtv_1.std()
print "Xtv_2.std() = ", Xtv_2.std()

max_diff = (abs(Xtv_1 - Xtv_2)).max()

print "Max deviance = ", max_diff

plotRegPath(f.add_subplot(413), Xtv_1 - Xtv_2)

a = f.add_subplot(414)
s1 = Xtv_1.reshape(Xtv_1.shape[0], -1).std(axis=1)
s2 = Xtv_2.reshape(Xtv_2.shape[0], -1).std(axis=1)

print "Ratio 1 = ", (s2 / s1).mean()
print "Ratio 2 = ", (s1 / s2).mean()

a.plot(range(Xtv_1.shape[0]), s1, label='xtv 1')
a.plot(range(Xtv_2.shape[0]), s2, label='xtv 2')
a.plot(range(Xtv_2.shape[0]), s2 / s1, label='ratio 1')
a.plot(range(Xtv_2.shape[0]), s1 / s2, label='ratio 2')
a.legend()

show()
