#!/use/bin/env python

import numpy as np
from pylatticeflow import calculate2dTV, calculate2dTVPath
from matplotlib.pylab import imread, figure, show
from matplotlib import collections
from itertools import product

def get_file(image_file, do_transform = True):
    Xo = imread(image_file)

    if not Xo.size:
        raise IOError("Error loading image %s." % image_file)

    X = (Xo.mean(axis=2) / Xo.max())

    if do_transform:
        X -= X.mean()
        X /= X.std()
    else:
        X /= X.max()

    return X

def get_ratio(X, lm):
    lambdas = np.linspace(0.05, lm, 10)

    Xtv_2 = calculate2dTVPath(X, lambdas)

    # Get the first round
    Xtv_1 = np.empty( (lambdas.size, X.shape[0], X.shape[1]) )

    for i, lm in enumerate(lambdas):
        Xtv_1[i,:,:] = calculate2dTV(X, lm)

    s1 = Xtv_1.reshape(Xtv_1.shape[0], -1).std(axis=1)
    s2 = Xtv_2.reshape(Xtv_2.shape[0], -1).std(axis=1)

    return (s1 / s2).mean();


x = np.linspace(0.1, 15, 50)

f = figure()
a = f.add_subplot(111)

X = get_file("benchmarks/images/sanity.png")
a.plot(x, [get_ratio(X, lm) for lm in x], '-b', label = "sanity")
a.plot(x, [np.sqrt(get_ratio(X, lm)) for lm in x], 'sr', label = "sanity-sqrt")

X = get_file("benchmarks/images/truffles-small.png")[::3,::3]
a.plot(x, [get_ratio(X, lm) for lm in x], ':g',  label = "truffles")
a.plot(x, [np.sqrt(get_ratio(X, lm)) for lm in x], 'dc', label = "truffles-sqrt")

a.set_xlim([0,x.max()])
a.legend()
show()




# image_file = "benchmarks/images/truffles.png"
#image_file = "benchmarks/images/truffles-small.png"
#image_file = "benchmarks/images/branches-small.png"
