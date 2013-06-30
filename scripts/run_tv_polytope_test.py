#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV
from matplotlib.pylab import imread, figure, show

image_file = "benchmarks/images/branches-small.png"

Xo = imread(image_file)

if not Xo.size:
    raise IOError("Error loading image %s." % image_file)

X = (Xo.mean(axis=2) / Xo.max())[:50,:50]

# X = np.linspace(0,1.25,6).reshape( (3,2) )

Xtv0 = calculate2dTV(X, 0)

assert abs(X - Xtv0).mean() <= 1e-4, abs(X - Xtv0).mean()

lambdas = np.exp(np.linspace(-5, 5, 1000))

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



f = figure()
a = f.add_subplot(221)
a.plot(Rtv[:,20, 20], Rtv[:,10, 20])
a = f.add_subplot(222)
a.plot(Rtv[:,20, 20], Rtv[:,20, 19])
a = f.add_subplot(223)
a.plot(Rtv[:,20, 20], Rtv[:,21, 20])
a = f.add_subplot(224)
a.plot(Rtv[:,20, 20], Rtv[:,20, 21])

show()




