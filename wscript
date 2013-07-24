#! /usr/bin/env python

# the following two variables are used by the target "waf dist"
VERSION='0.0.1'
APPNAME='latticeqpb'

top = '.'
out = 'build'

import waflib.Configure
from glob import glob
import os
from os.path import join, abspath
#waflib.Configure.autoconfig = True
import install as W
import install.code_setup as K
import imp
import sys
from ctypes import cdll

__solver_list = ["LatticeEnergyMinimizer",
                 "LatticeLevelReductions"
                 ]

################################################################################
# Basic configuration information; detects possible flags

def options(ctx):
    ctx.load('compiler_cxx')

    ctx.add_option('-g', '--debug', dest='debug', action='store_true', default=False,
                   help='Turn on debug mode with optimization and numerous checks.')

    ctx.add_option('--debug-full', dest='debug_full', action='store_true', default=False,
                   help='Turn on full debug mode with no optimization and numerous checks.')

    ctx.add_option('--enable-pr-checks', dest='enable_pr_checks', action='store_true', default=False,
                   help='Enable checks for the push-relable engine.')

    ctx.add_option('', '--fast', dest='fast', action='store_true', default=False,
                   help='Turn on full optimization options at the expense of compile time.')

    ctx.add_option('', '--quick-config', dest='quick_config', action='store_true', default=False,
                   help='Assume the important options just work and skip tests.')

    ctx.add_option('--no-inline', dest='no_inline', action='store_true', default=False,
                   help='Prevent inlining of functions for profiling purposes.')

    ctx.add_option('--1d', dest='enable_1d', action='store_true', default=False,
                   help='Enable generation of 1d kernels')
    ctx.add_option('--max-1d-radius', dest='max_1d_radius', default=8,
                   help='The maximum radius to use in creating 1d edge kernels; default = 4.')

    ctx.add_option('--no-2d', dest='enable_2d', action='store_false', default=True,
                   help='Disable generation of 2d kernels')
    ctx.add_option('--max-2d-radius', dest='max_2d_radius', default=4,
                   help='The maximum radius to use in creating 2d edge kernels; default = 4.')

    ctx.add_option('--3d', dest='enable_3d', action='store_true', default=False,
                   help='Enable generation of 3d kernels')
    ctx.add_option('--max-3d-radius', dest='max_3d_radius', default=2,
                   help='The maximum radius to use in creating 3d edge kernels; default = 2.')

    ctx.add_option('--4d', dest='enable_4d', action='store_true', default=False,
                   help='Enable generation of 4d kernels')
    ctx.add_option('--max-4d-radius', dest='max_4d_radius', default=1,
                   help='The maximum radius to use in creating 4d edge kernels; default = 1.')

    ctx.add_option('--all-kernels', dest='all_kernels', action='store_true', default=False,
                       help='Enable creation of 1d, 2d, 3d, and 4d kernels.')    

    ctx.add_option('--test', dest='test_mode', action='store_true', default=False,
                       help= ('Test/debugging mode: all debug flags, '
                              'plus only one kernel for debugging purposes '
                              '(can be overridden by fast).'))

    ctx.add_option('--enable-bkgc', dest='enable_bkgc', action='store_true', default=False,
                   help = ('Enable checks and comparison with the Boykov-Kolmogorov '
                           'graphcuts code for lattice.'))

    ctx.add_option('--no-python', dest='python', action='store_false', default=True,
                   help='Disable compilation of the python extension module.')

    ctx.add_option('--benchmark', dest='benchmark', default="quick",
                   help='Specify the type of benchmark to run (default: quick).')

    ctx.add_option('--benchmark-checks', dest='benchmark_check', action='store_true',
                   default=False,
                   help='Specify whether to check the output solutions of the benchmarks.')


    W.cython.options(ctx)

def __process_cascading_options(ctx):

    if ctx.options.all_kernels:
        ctx.options.enable_1d = True
        ctx.options.enable_2d = True
        ctx.options.enable_3d = True
        ctx.options.enable_4d = True

    if ctx.options.test_mode:
        ctx.options.enable_1d = False
        ctx.options.enable_2d = True
        ctx.options.enable_3d = False
        ctx.options.enable_4d = False

        ctx.options.max_2d_radius = 1
        ctx.options.enable_bkgc = True

        ctx.options.benchmark = "test"

    # Configure these options
    if not ctx.options.enable_1d:
        ctx.options.max_1d_radius = 0
    if not ctx.options.enable_2d:
        ctx.options.max_2d_radius = 0
    if not ctx.options.enable_3d:
        ctx.options.max_3d_radius = 0
    if not ctx.options.enable_4d:
        ctx.options.max_4d_radius = 0

    

def configure(ctx):

    # Prefilter the options
    __process_cascading_options(ctx)
    
    ctx.load('compiler_cxx')
    
    W.setupCompilerFlags(ctx, cxx11support = True)
    
    W.checkCXX11features(ctx, 
                         ['auto',
                          'lamda',
                          'decltype',
                          'nullptr',
                          'rvalue reference',
                          "static_assert",
                          'range based for loop',
                          'initializer list',
                          'variadic templates',
                          'SFINAE',
                          'variadic macro support',
                          'constexpr',
                          'template arguments',
                          'template aliasing'
                          ])

    # Set up python if need be
    if ctx.options.python:
        ctx.load('python')
        ctx.check_python_headers()

        ctx.load("cython")
        W.cython.configure(ctx)

    ctx.env.mode = "release"

    if ctx.options.debug:
        ctx.env.mode = "debug"
    if ctx.options.debug_full:
        ctx.env.mode = "debug-full"
    if ctx.options.fast:
        ctx.env.mode = "fast"

    cxxflags  = ctx.env.cxxflags[ctx.env.mode]
    linkflags = ctx.env.linkflags[ctx.env.mode]

    if ctx.options.no_inline:
        cxx_flags.append("-fno-inline")

    if ctx.options.enable_bkgc:
        cxxflags.append("-DENABLE_BOYKOV_KOLMOGOROV_GC_CODE")
        linkflags.append("-DENABLE_BOYKOV_KOLMOGOROV_GC_CODE")

    if ctx.options.enable_pr_checks:
        cxxflags.append("-DENABLE_PR_CHECKS")
        linkflags.append("-DENABLE_PR_CHECKS")

    ctx.env.append_value('CXXFLAGS', cxxflags)
    ctx.env.append_value('LINKFLAGS', linkflags)

    ################################################################################
    # Set up the solver list
    ctx.env.solver_list = __solver_list
    
    if ctx.options.enable_bkgc:
        ctx.env.solver_list.append("BKGCEnergyMinimizer")

    ################################################################################
    # Setting up the stuff with the kernel
        
    ctx.env.enable_1d = ctx.options.enable_1d
    ctx.env.enable_2d = ctx.options.enable_2d
    ctx.env.enable_3d = ctx.options.enable_3d
    ctx.env.enable_4d = ctx.options.enable_4d

    ctx.env.max_1d_radius = ctx.options.max_1d_radius
    ctx.env.max_2d_radius = ctx.options.max_2d_radius
    ctx.env.max_3d_radius = ctx.options.max_3d_radius
    ctx.env.max_4d_radius = ctx.options.max_4d_radius

    # First, set up the kernel files
    K.generateKernelSources(ctx)

    # Now add all of the directories under solver to the path
    for d in set(abspath(d) for d, dl, fl in os.walk('solver')):
        ctx.env.append_value('INCLUDES', d)

def build(ctx):

    ctx.program(source       = 'solver/benchmark.cpp', 
                target       = 'latticeqbp_random_benchmark',
                includes     = ['solver'], 
                lib          = ['m'])

    main_source_files = glob(join(K.kernel_explicit_instantiation_directory, '*.cpp'))
    main_source_files.append(K.factory_write_file)

    ctx.stlib(source=main_source_files,
              target='stlatticeflow',
              lib          = ['m'],
              install_path = join(ctx.options.prefix, 'lib'))

    ctx.shlib(source=[],
              target='latticeflow',
              install_path = join(ctx.options.prefix, 'lib'),
              lib          = ['m'],
              use = ['stlatticeflow']
              )

    if ctx.options.python:
        ctx(features = 'cxx cxxshlib pyext cython',
            source   = ['extensions/python/pylatticeflow.pyx'],
            target   = 'pylatticeflow',
            use      = ['stlatticeflow']
            )

def benchmark(ctx):
    
    __process_cascading_options(ctx)

    if ctx.options.benchmark is not None:
        ctx.env.benchmark = ctx.options.benchmark

    if ctx.options.benchmark_check:
        ctx.env.benchmark_check = ctx.options.benchmark_check

    try:
        py_mod_info = imp.find_module('pylatticeflow', [abspath(out)])
    except ImportError:
        ctx.fatal("Not finding compiled pylatticeflow module "
                  "in build directory; has build been run?")

    # Find the liblatticeflow library
    pylatticeflow = imp.load_module('pylatticeflow', *py_mod_info)

    def runBasicBenchmark(solvers, image_files, kernels):

        print_list = []

        image_names = [im.replace('benchmarks/images/', '') for im in image_files]

        kernel_padding = 2 + max(len("{Kernel:%s}" % k) for k in kernels)
        solver_padding = 2 + max(len("{Solver:%s}" % s) for s in solvers)
        image_padding = 2 + max(len("{Image:%s}" % im) for im in image_names)

        for k in kernels:

            ks = "{Kernel:%s}" % k
            ks += " "*(kernel_padding - len(ks))

            for s in solvers:

                ss = "{Solver:%s}" % s
                ss += " "*(solver_padding - len(ss))

                total_time = 0

                for im_name, im_file in zip(image_names, image_files):

                    try:
                        time = pylatticeflow._runBenchmark(im_file, k, s, ctx.env.benchmark_check)
                    except ValueError, ve:
                        print "Error running test: %s" % str(ve)
                        continue

                    total_time += time

                    im_string = "{Image:%s}" % im_name
                    im_string += " "*(image_padding - len(im_string))

                    ps = "%s %s %s {Time:%f}" % (ks, ss, im_string, time)
                    print ps

                    print_list.append(ps)

                if len(image_names) > 1:
                    ps = "%s %s {Image:Average} {Time:%f}" % (ks, ss, total_time / len(image_names))
                    print ps
                    print_list.append(ps + "\n")

        print "#"*45
        print "\n".join(print_list)


    ################################################################################

    solver_list = (["LatticeEnergyMinimizer"]
                   + (["BKGCEnergyMinimizer"] if ctx.env.enable_bkgc else []))

    if ctx.env.benchmark == "quick":

        runBasicBenchmark(
            solvers = solver_list,
            image_files = ["benchmarks/images/branches.png"],
            kernels = ["Full2d_12", "Full2d_4"]
            )

    if ctx.env.benchmark == "full":

        runBasicBenchmark(
            solvers = solver_list,
            image_files = glob("benchmarks/images/*.png"),
            kernels = ["Full2d_12", "Full2d_4", "Full2d_48"],
            )

    elif ctx.env.benchmark == "test":

        runBasicBenchmark(
            solvers = solver_list,
            image_files = ["benchmarks/images/branches-small.png"],
            kernels = ["Full2d_4"])

    elif ctx.env.benchmark == "sanity":

        runBasicBenchmark(
            solvers = ["LatticeEnergyMinimizer"],
            image_files = ["benchmarks/images/sanity.png"],
            kernels = ["Full2d_4"])
    else:

        ctx.fatal("Benchmark %s not recognized; available benchmarks are: test, quick, quick3d, full"
                  % ctx.options.benchmark)


from waflib.Build import BuildContext
class benchmark_(BuildContext):
    cmd = 'benchmark'
    fun = 'benchmark'
