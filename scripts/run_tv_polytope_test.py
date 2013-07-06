#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV
from matplotlib.pylab import imread, figure, show
from matplotlib import collections
from itertools import product

plot_type = "levels2"

image_file = "benchmarks/images/truffles-small.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max()) #[::2, ::2]

X -= X.mean()

# X = np.linspace(0,1.25,6).reshape( (3,2) )

# Xtv0 = calculate2dTV(X, 0)

# assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

lambdas = np.linspace(0, .1, 100)

Xtv = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )
Rtv = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )

slmdas = lambdas

for i, lm in enumerate(lambdas):
    Xtv[i,:,:] = calculate2dTV(X, lm)
    # Rtv[i,:,:] = (Xtv[i,:,:] - X) / (1e-32 + lambdas[i])

print "Done calculating models."


# import pandas as pd
# import matplotlib.pyplot as plt 

# df = pd.DataFrame(Rtv[:,25:50:5,:25:30].reshape(Rtv.shape[0], -1))
# axes = pd.tools.plotting.scatter_matrix(df)
# plt.tight_layout()
# plt.show()

if plot_type == "pairs":

    def plot_rtv(a, xt, yt):

        x, y = 10000*Rtv[:,xt[0], xt[1]], 10000*Rtv[:,yt[0], yt[1]]

        min_x, max_x = x.min(), x.max()
        min_y, max_y = y.min(), y.max()

        cx = (max_x + min_x) / 2
        cy = (max_y + min_y) / 2

        w = max(max_x - min_x, max_y - min_y, 2*abs(cx), 2*abs(cy))

        M = max(abs(min_x), abs(max_x), abs(min_y), abs(max_y))

        a.plot(x, y, '+-b')
        a.plot([-M, M], [-M, M], '--k')
        a.plot([x[0]], [y[0]], 'xr', markersize = 30)
        a.set_xlim([cx - w/2, cx + w/2])
        a.set_ylim([cy - w/2, cy + w/2])
        a.set_title(str(xt) + ":" + str(yt)
                    + (": v=%0.5f,%0.5f" % (X[xt[0], xt[1]], X[yt[0], yt[1]])))


    fv = [
        ( (30, 30), (30, 31) ),
        ( (30, 30), (30, 32) ),
        ( (30, 30), (30, 33) ),
        ( (30, 30), (30, 34) ),
        # ( (30, 34), (30, 35) ),
        # ( (30, 35), (30, 36) ),
        # ( (30, 36), (30, 37) ),
        # ( (30, 37), (30, 38) ),
        # ( (31, 38), (31, 39) ),
        ( (30, 30), (31, 31) ),
        ( (30, 30), (31, 32) ),
        ( (30, 30), (31, 33) ),
        ( (30, 30), (31, 34) ),
        # ( (31, 34), (31, 35) ),
        # ( (31, 35), (31, 36) ),
        # ( (31, 36), (31, 37) ),
        # ( (31, 37), (31, 38) ),
        # ( (31, 38), (31, 39) ),
        ( (30, 30), (32, 31) ),
        ( (30, 30), (32, 32) ),
        ( (30, 30), (32, 33) ),
        ( (30, 30), (32, 34) ),
        # ( (32, 34), (32, 35) ),
        # ( (32, 35), (32, 36) ),
        # ( (32, 36), (32, 37) )
        ]


    f = figure()

    for i, ft in enumerate(fv):
        a = f.add_subplot(4,3,i + 1)
        plot_rtv(a, *ft)

    show()


if plot_type == "single":

    xi, yi = 5, 5
    width = 5

    f = figure()
    a = f.add_subplot(111)

    # x -= 1000*Rtv.reshape(Rtv.shape[0], -1).mean(axis=1)

    # for i, j in product(range(xi-width, xi + width),
    #                     range(yi-width, yi + width) ):

    for xi, yi in product(range(1, Rtv.shape[1]-1), range(1, Rtv.shape[2] -1)):

        for i, j in [(xi + dx, yi + dy)
                     for dx, dy in [(-1, 0), (1,0), (0,-1), (0,1)]]:

            x = 1000*Rtv[:,xi, yi]

            y = 1000*Rtv[:,i, j]# - x

            # a.plot([x[0]], [y[0]], 'xr', markersize = 30)
            y = 1000*Rtv[:,i, j] - x

            if y.min() < -1 and y.max() > 1:
                a.plot(y[::-1], label="%d,%d -> %d,%d" % (xi, yi, i,j) )

    a.plot([0,0], [len(x), 0], '--k')
    a.legend()
    show()

if plot_type == "levels":

    f = figure()
    a = f.add_subplot(111)

    #for xi, yi in product(range(0, Rtv.shape[1], 5), range(0, Rtv.shape[2], 5)):
    for xi, yi in product(range(0, 10), range(0, 10)):
        
        y = 1000*Rtv[:,xi, yi]

        a.plot(lambdas[::-1], y[::-1], label="%d,%d" % (xi, yi))

    a.legend()
    show()


if plot_type == "levels2":

    f = figure()
    a = f.add_subplot(111)

    xc, yc = 50, 50
    lmi = int(lambdas.size * 0.25)

    col = []

    ymin, ymax = 0,0

    for xi, yi in product(range(0, Rtv.shape[1]), range(0, Rtv.shape[2])):
    #for xi, yi in product(range(0, 10), range(0, 10)):
        
        y = lambdas * Xtv[:,xi, yi]

        yref = lambdas * Xtv[:, xc, yc]

        if abs(y[:lmi] - yref[:lmi]).mean() < 0.00001:
            col.append( np.vstack([lambdas, y]).T )
            ymin = min(ymin, y.min())
            ymax = max(ymax, y.max())

    a.add_collection(collections.LineCollection(col))
    a.legend()
    a.set_xlim([0,lambdas.max()])
    a.set_ylim([ymin, ymax])
    show()
