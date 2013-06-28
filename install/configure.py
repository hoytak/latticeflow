################################################################################
# Configuration stuff for the compilers

from copy import copy

def addOption(ctx, *options):

    for option in options:

        if type(option) is str:

            ctx.env.append_value('CXXFLAGS', option)

        elif type(option) in [list, tuple]:

            for opt in option:
                ctx.env.append_value('CXXFLAGS', opt)

        else:
            raise TypeError("Options specification has the wrong type.")

def compilerOptionWorks(ctx, option, base_options = [], must_work = False):
    ctx.env.stash()

    # ctx.start_msg("Checking %s option '%s': "
    #               % ("required" if must_work else "optional", str(option)))

    addOption(ctx, *base_options)
    addOption(ctx, '-Werror', option)

    try:
        ctx.check_cxx(
            msg = ("Checking %s support for compiler option %s"
                   % ("required" if must_work else "optional", option)),
            define_name = str(option),
            fragment    = "int main() {return 0; }",
            execute     = True)

        return True

    except Exception:
        if must_work:
            ctx.fatal("Required compiler option %s doesn't work \n  > Base options = %s."
                      % (option, ','.join(str(opt) for opt in base_options)))
        return False

    finally:
        ctx.env.revert()

def filterOptions(ctx, options, base_options = [], all_must_work = False):

    base_options = copy(base_options)

    ok_options = []

    for option in options:
        if compilerOptionWorks(ctx, option,
                               base_options = base_options,
                               must_work = all_must_work):

            base_options.append(option)
            ok_options.append(option)

    return ok_options


def chooseOption(ctx, *options):
    for opt in options:
        if compilerOptionWorks(ctx, opt):
            return opt

    ctx.fatal("One of the options %s needs to be present." % (', '.join(options)) )
