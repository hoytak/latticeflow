#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV
from matplotlib.pylab import imread, figure, show

image_file = "benchmarks/images/branches-small.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())[::2, ::2]

# X = np.linspace(0,1.25,6).reshape( (3,2) )

Xtv0 = calculate2dTV(X, 0)

assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

lambdas = np.exp(np.linspace(-5, 16, 500))
lambdas[0] = 1e-8

Xtv = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )
Rtv = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )

for i, lm in enumerate(lambdas):
    Xtv[i,:,:] = calculate2dTV(X, lm)
    Rtv[i,:,:] = (Xtv[i,:,:] - X) / lambdas[i]


# import pandas as pd
# import matplotlib.pyplot as plt 

# df = pd.DataFrame(Rtv[:,25:50:5,:25:30].reshape(Rtv.shape[0], -1))
# axes = pd.tools.plotting.scatter_matrix(df)
# plt.tight_layout()
# plt.show()

def plot_rtv(a, xt, yt):

    x, y = Rtv[:,xt[0], xt[1]], Rtv[:,yt[0], yt[1]]

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
    a.set_title(str(xt) + ":" + str(yt))


fv = [
    ( (30, 30), (30, 31) ),
    ( (30, 31), (30, 32) ),
    ( (30, 32), (30, 33) ),
    ( (30, 33), (30, 34) ),
    ( (30, 34), (30, 35) ),
    ( (30, 35), (30, 36) ),
    ( (30, 36), (30, 37) ),
    ( (30, 37), (30, 38) ),
    ( (30, 38), (30, 39) )
    ]
      

f = figure()

for i, ft in enumerate(fv):
    a = f.add_subplot(3,3,i + 1)
    plot_rtv(a, *ft)

show()




