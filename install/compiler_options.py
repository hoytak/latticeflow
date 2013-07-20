from configure import chooseOption, filterOptions, addOption
import cpp11x
import os

def setupCompilerFlags(ctx, cxx11support = False):

    compiler = "g++"

    ############################################################
    #
    if cxx11support:
        cpp11x.setupCXX11(ctx)

    ############################################################
    # Set up the compiler options to include libraries

    include_path = os.getenv("INCLUDE_PATH")

    if include_path is not None:
        ip_includes = [s.strip() for s in include_path.split(':') if s.strip()]
        ctx.env.append_value('INCLUDES', ip_includes)

    ############################################################
    #  Set up the options for gcc
    if compiler in ["g++", "gcc", "clang", "clang++"]:

        compiler_common_options = [cpp11x.cxx11_flag, '-Wall', '-g', '-fpic',
                                   '-Wno-unused-local-typedefs']
        linker_common_options   = [cpp11x.cxx11_flag, '-Wall', '-g', '-fpic']

        compiler_full_debug_options = ['-O0', '-Wno-strict-aliasing', '-fno-inline']
        linker_full_debug_options = []

        compiler_debug_options = [
            '-O2', '-Wno-strict-aliasing', '-UNDEBUG']
        linker_debug_options = ['-UNDEBUG']

        compiler_release_options = [
            '-O3', '-funroll-all-loops', '-ffast-math', '-DNDEBUG']
        linker_release_options  = ['-O3', '-DNDEBUG']

        compiler_super_options = [
            '-march=native',
            '-fipa-struct-reorg',
            '-fipa-cp',
            '-ftree-loop-linear',
            '-mfpmath=sse',
            '-ffast-math',
            '-fstrict-aliasing',
            '-floop-interchange',
            '-floop-block',
            '-floop-strip-mine',
            '-fno-common',
            '-fmodulo-sched-allow-regmoves',
            '-fmodulo-sched',
            '-finline-limit=20000',
            '-fgraphite-identity',
            '-ftree-loop-distribution',
            '-fivopts',
            '-ftracer',
            '-fvariable-expansion-in-unroller',
            '-freorder-blocks-and-partition',
            '-freciprocal-math',
            '-fprefetch-loop-arrays',
            '-funroll-loops',
            '-fweb',
            '-funsafe-math-optimizations',
            '-ffinite-math-only',
            '-fearly-inlining',
            '-fpeel-loops',
            '-fbranch-target-load-optimize',
            '-fgcse-sm',
            '-fgcse-las',
            '-funsafe-loop-optimizations',
            ('--param', 'inline-unit-growth=400'),
            ('--param', 'max-inline-insns-auto=1000'),
            ('--param', 'ipcp-unit-growth=1000'),
            ('--param', 'max-inline-recursive-depth=32'),
            ('--param', 'min-insn-to-prefetch-ratio=1'),
            ('--param', 'prefetch-min-insn-to-mem-ratio=1'),
            ('--param', 'simultaneous-prefetches=32')
            ]

        linker_super_options = []
    else:
        ctx.fatal("Compiler '%s' not recognized." % compiler)

    ################################################################################
    # Now will need to do the same for VC++

    # TODO

    ################################################################################
    # Now check each of the options

    # filterOptions(ctx, compiler_common_options, all_must_work = True)
    # filterOptions(ctx, linker_common_options, all_must_work = True)

    def compiler_fltr(args):
        return filterOptions(ctx, args,
                             base_options=compiler_common_options,
                             all_must_work = False)

    def linker_fltr(args):
        return filterOptions(ctx, args,
                             base_options=linker_common_options,
                             all_must_work = False)

    compiler_debug_options = \
        compiler_fltr(compiler_debug_options)
    linker_debug_options = \
        linker_fltr(linker_debug_options)

    compiler_full_debug_options = \
        compiler_fltr(compiler_full_debug_options)
    linker_full_debug_options = \
        linker_fltr(linker_full_debug_options)

    compiler_release_options = \
        compiler_fltr(compiler_release_options)
    linker_release_options = \
        linker_fltr(linker_release_options)

    compiler_super_options = \
        compiler_fltr(compiler_super_options)
    linker_super_options = \
        linker_fltr(linker_super_options)

    ctx.env.cxxflags = {}
    ctx.env.linkflags = {}

    ctx.env.cxxflags["debug"]  = compiler_common_options + compiler_debug_options
    ctx.env.linkflags["debug"] = linker_common_options   + linker_debug_options

    ctx.env.cxxflags["debug-full"]  = ctx.env.cxxflags["debug"]  + compiler_full_debug_options
    ctx.env.linkflags["debug-full"] = ctx.env.linkflags["debug"] + linker_full_debug_options

    ctx.env.cxxflags["release"]  = compiler_common_options + compiler_release_options
    ctx.env.linkflags["release"] = linker_common_options   + linker_release_options

    ctx.env.cxxflags["fast"]  = ctx.env.cxxflags["release"]  + compiler_super_options
    ctx.env.linkflags["fast"] = ctx.env.linkflags["release"] + linker_super_options

    def optAsString(opt):
        return opt if type(opt) is str else " ".join(opt)

    ctx.msg("Compiler flags for debug mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.cxxflags["debug"]))
    ctx.msg("Link flags for debug mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.linkflags["debug"]))

    ctx.msg("Compiler flags for full debug mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.cxxflags["debug-full"]))
    ctx.msg("Link flags for full debug mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.linkflags["debug-full"]))

    ctx.msg("Compiler flags for release mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.cxxflags["release"]))
    ctx.msg("Link flags for release mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.linkflags["release"]))

    ctx.msg("Compiler flags for fast mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.cxxflags["fast"]))
    ctx.msg("Link flags for fast mode",
            "\n" + " ".join(optAsString(s) for s in ctx.env.linkflags["fast"]))
