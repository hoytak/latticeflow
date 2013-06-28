from pylatticeflow import LatticeFlow, _runBenchmark
from scipy.ndimage import imread



def runBenchmark(image_file, dimensions, kernel, solver, test_cut):

    d = {}
