import math
import re
from itertools import product
from collections import defaultdict
from os.path import exists, join, split, abspath
from os import makedirs, getcwd
from hashlib import md5
from math import floor, ceil, sqrt
from shutil import rmtree

# Min included kernel templates
min_kernel_dimensions = {1: 0, 2:4, 3:2, 4:0}

# These variables the paths to the rest of the project.
kernel_write_file = 'solver/kernels/kernels.hpp'
kernel_explicit_instantiation_directory = 'solver/kernels/explicit_instantiations/'
instantiation_suppression_write_file = 'solver/kernels/instantiation_suppressions.hpp'
factory_write_file = 'solver/interfacing/factory.cpp'

documentation_directory = 'doc/kernels/'
documentation_file_name = '%(dimension)ddkernels.rst'

################################################################################

graph_scale = 0.4

################################################################################

def load_template(s):
    return open(join("install/templates", s), "r").read()

class Info: pass

T = Info()

T.header = load_template("header.txt")

############################################################
# Names and small specifics
T.kernel_name = "%(kernel_type)s%(dimension)dd_%(n_edges)s"
T.kernel_header_name = T.kernel_name

############################################################
# Files and more complicated templates

T.kernel_file                     = load_template("kernels/file.hpp")
T.kernel_code                     = load_template("kernels/code.hpp")
T.kernel_duplicate_code           = load_template("kernels/duplicate_code.hpp")

T.kernel_documentation_file       = load_template("kernels/documentation_file.rst")
T.kernel_documentation_header     = load_template("kernels/documentation_header.rst")
T.kernel_documentation            = load_template("kernels/documentation.rst")
T.kernel_duplicate_documentation  = load_template("kernels/duplicate_documentation.rst")

T.graphviz_node                   = load_template("graphviz/node.dot")
T.graphviz_edge                   = load_template("graphviz/edge.dot")
T.graphviz_grid_edge              = load_template("graphviz/grid_edge.dot")

T.explicit_instantiation_file     = load_template("explicit_instantiation_file.cpp")
T.explicit_instantiation_code     = load_template("explicit_instantiation_code.cpp")
T.instantiation_suppression_code  = load_template("factory/instantiation_suppression_code.hpp")

T.factory_file                    = load_template("factory/file.cpp")
T.factory_creation_function       = load_template("factory/creation_function.cpp")
T.factory_kernel_dimension_entry  = load_template("factory/kernel_dimension.cpp")
T.factory_map_entry               = load_template("factory/map_entry.cpp")

################################################################################
# The code that controls which kernels get generated

def generateKernelSources(ctx):

    def getRadiusIntervals(dimension, max_radius):
        return sorted(d for d in set( sqrt(sum(x**2 for x in loc))
                                      for loc in product(*([range(0, max_radius+1)] * dimension)))
                      if 0 < d <= max_radius + 1e-8)

    kf = SourceFileSetup(ctx)

    ############################################################
    # 1 Dimension
    if ctx.env.enable_1d or min_kernel_dimensions[1] != 0:

        ctx.start_msg("Preparing 1d kernels... ")

        for radius in xrange(1, max(ctx.env.max_1d_radius, min_kernel_dimensions[1])):
            kf.addKernel("Full", 1, radius,
                         "Connection to all nodes within %s of each other." % radius,
                         radius > ctx.env.max_1d_radius
                         )

        ctx.end_msg("Done.")

    ############################################################
    # 2 Dimensions

    if ctx.env.enable_2d or min_kernel_dimensions[2] != 0:
        ctx.start_msg("Preparing 2d kernels... ")

        for radius in getRadiusIntervals(2, max(ctx.env.max_2d_radius, min_kernel_dimensions[2])):

            kf.addKernel("Full", 2, radius,
                         "Connection to all nodes within %0.2f of each other." % radius,
                         radius > ctx.env.max_2d_radius
                         )

            kf.addKernel("Star", 2, radius,
                         "Connection to all nodes within %0.2f of each other not overlapping shorter paths." % radius,
                         radius > ctx.env.max_2d_radius
                         )

        ctx.end_msg("Done.")

    ############################################################
    # 3 Dimensions

    if ctx.env.enable_3d or min_kernel_dimensions[3] != 0:

        ctx.start_msg("Preparing 3d kernels... ")

        for radius in getRadiusIntervals(3, max(ctx.env.max_3d_radius, min_kernel_dimensions[3])):
            kf.addKernel("Full", 3, radius,
                         "Connection to all nodes within %0.2f of each other." % radius,
                         radius > ctx.env.max_3d_radius
                         )

            kf.addKernel("Star", 3, radius,
                         "Connection to all nodes within %0.2f of each other not overlapping shorter paths." % radius,
                         radius > ctx.env.max_3d_radius
                         )

        ctx.end_msg("Done.")

    ############################################################
    # 4 Dimensions

    if ctx.env.enable_4d or min_kernel_dimensions[4] != 0:

        ctx.start_msg("Preparing 3d kernels... ")

        for radius in getRadiusIntervals(4, max(ctx.env.max_4d_radius, min_kernel_dimensions[4])):
            kf.addKernel("Full", 4, radius,
                         "Connection to all nodes within %0.2f of each other." % radius,
                         radius > ctx.env.max_4d_radius
                         )


            kf.addKernel("Star", 4, radius,
                         "Connection to all nodes within %0.2f of each other not overlapping shorter paths." % radius,
                         radius > ctx.env.max_4d_radius
                         )

        ctx.end_msg("Done.")

    kf.write()


################################################################################
# Now the code to generate certain kernels

class KernelEdgeGenerator:

    def __getGeocutEdgeWeights(self, dimension, edges):

        if dimension == 2 and len(edges) >= 4:

            angles = [math.atan2(e[1], e[0]) % (2*math.pi) for e in edges]

            # Verify that it's ok; if not, then just say so

            if len(set(angles)) != len(edges):
                return (False, [1]*len(edges))
            
            edge_info = sorted(zip(angles, edges))
            angles = [aw for aw, e in edge_info]

            def get_aw(idx):
                n = len(angles)
                a1 = angles[ (idx - 1) % n]
                a2 = angles[ (idx + 1) % n]

                if a2 < a1:
                    a2 += 2*math.pi

                return (a2 - a1) / 2

            angle_widths = [get_aw(i) for i in xrange(len(angles))]

            assert abs(sum(angle_widths) - 2*math.pi) < 1e-8, sum(angle_widths)

            weight_map = dict(
                (e, aw / (2.0 * math.sqrt(e[0]*e[0] + e[1]*e[1])))
                for aw, (ang, e) in zip(angle_widths, edge_info))

            return (True, [weight_map[e] for e in edges])
        else:
            return (False, [1]*len(edges))
        

    def Full(self, dimension, radius):
        """
        Makes things in a solid sphere pattern.
        """

        d, h = dimension, radius

        lr = []

        def is_valid(x):
            return

        iters = [range(int(floor(-h)), int(ceil(h))+1)] * d

        for x in product( *iters):
            if 0 < sum(xe**2 for xe in x) <= radius**2:
                lr.append(tuple(x))

        lr = sorted(lr)

        R = Info()

        R.edges = lr
        R.is_geocut_applicable, R.geocut_edge_weights = \
            self.__getGeocutEdgeWeights(dimension, lr)

        return R

    def Star(self, dimension, radius):
        """
        Makes things in a sphere pattern, but with no lines crossing
        other lines.
        """

        lr = set(self.Full(dimension,radius).edges)

        # print "$$$$$$$$$$$$$$$$$$ STAR"

        lr_s = sorted(lr, key = lambda edge: sum(u*u for u in edge))
        angles = [tuple(math.atan2(u, e[0]) for u in e[1:]) for e in lr_s]

        # print "\n".join(str(edge + (ang,)) for edge, ang in zip(lr_s, angles))

        seen_angles = set()
        del_set = set()
        
        for ang, e in zip(angles, lr_s):
            if ang in seen_angles:
                lr.remove(e)
            else:
                seen_angles.add(ang)

        edges = sorted(lr)

        # print "geocut!!!!!!!!!!!!!!!!"
        
        R = Info()
        R.edges = edges
        R.is_geocut_applicable, R.geocut_edge_weights = \
            self.__getGeocutEdgeWeights(dimension, lr)

        # print "<<<<<<<<<<<<<<<<<<<<<< DONE" 

        return R

################################################################################
# Now the meet of the templating

class SourceFileSetup:

    def __init__(self, ctx):

        base_dir = getcwd()

        self.kernel_write_file = abspath(join(base_dir, kernel_write_file))
        self.documentation_directory = abspath(join(base_dir, documentation_directory))
        self.kernel_explicit_instantiation_directory = \
            abspath(join(base_dir, kernel_explicit_instantiation_directory))
        self.instantiation_suppression_write_file = \
            abspath(join(base_dir, instantiation_suppression_write_file))
        self.factory_write_file = abspath(join(base_dir, factory_write_file))

        assert exists(split(self.kernel_write_file)[0]), \
               ("Kernal write directory '%s' does not exist; running script from base of project tree?" % split(self.kernel_write_file)[0])

        assert exists(self.documentation_directory), \
               ("Kernal documentation directory '%s' does not exist; running script from base of project tree?" % self.documentation_directory)

        if exists(self.kernel_explicit_instantiation_directory):
            rmtree(self.kernel_explicit_instantiation_directory)
        makedirs(self.kernel_explicit_instantiation_directory)

        self.solver_list = ctx.env.solver_list

        self.kernel_source_list = []
        self.instantiation_file_list = []

        self.documentation_file_list = defaultdict(lambda: [])
        self.documentation_file_list.update(dict( ((i, documentation_file_name % {"dimension" : i}), []) for i in [1,2,3,4]))

        self.factory_generator_map = []
        self.kernel_list = []
        self.factory_creation_functions = []
        self.factory_kernel_dimension_list = []
        self.instantiation_suppression_list = []

        self.seen_kernels = {}
        self.seen_kernel_names = set()

        self.kernel_dimensions_considered = set()

        self.ctx = ctx

    def getDictionary(self, ktype, dimension, radius, description,
                      template_only, **kwargs):

        keg = KernelEdgeGenerator()

        edge_info = getattr(keg, ktype)(dimension, radius, **kwargs)

        edge_list = edge_info.edges

        if len(edge_list) == 0:
            return True, None

        ########################################
        # fill in the basic info

        d = {}

        d['dimension'] = dimension
        self.kernel_dimensions_considered.add(dimension)

        d['radius'] = radius

        main_r = int(floor(ceil(100.0*radius) / 100.0))
        remainder = int(ceil( (100.0*radius))) % 100

        d['radius_string'] = ('%dp%02d' % (main_r, remainder)
                              if remainder != 0 else
                              '%d' % main_r)
        ################################################################################
        # Create the edge string

        edge_list_str = '{'

        for i, x in enumerate(edge_list):
            edge_list_str += ('{' + ','.join(str(xe) for xe in x)
                              + ('},' if i != len(edge_list) - 1 else '}'))

            if (i+1) % (dimension * 2) == 0:
                edge_list_str += '\n '

        edge_list_str += '}'


        d['edge_list'] = edge_list_str
        d['is_geocut_applicable'] = 1 if edge_info.is_geocut_applicable else 0
        d['geocut_edge_weights'] = (
            '{' + ', '.join( ('%1.16f' % w if w != 1 else "1.0")
                            for w in edge_info.geocut_edge_weights) + '}')
        
        d['n_edges'] = len(edge_list)

        ############################################################
        #

        d["kernel_type"] = ktype

        if "name" in kwargs:
            d['kernel_name'] = name = kwargs["name"]

            if name in self.seen_kernel_names:
                self.ctx.fatal("Configuration Error: Given kernel name '%s' already in use!!" % name)

            if re.match("^[a-zA-Z_][a-zA-Z0-9_]*$", name) == None:
                self.ctx.fatal("Configuration Error: Given kernel name '%s' not a vald C++ identifier." % name)

        else:
            d['kernel_name'] = T.kernel_name % d

        d['kernel_header_name'] = T.kernel_header_name % d
        d['documentation_file'] = documentation_file_name % d

        kernel_key = frozenset(edge_list)

        if kernel_key in self.seen_kernels:
            if d['kernel_name'] in self.seen_kernel_names:
                return False, None

            d_base = self.seen_kernels[kernel_key].copy()
            d["reference_kernel_name"] = d_base["kernel_name"]
            d["reference_kernel_type"] = d_base["kernel_type"]
            d_base.update(d)
            d = d_base

            d["documentation"] = T.kernel_duplicate_documentation % d

            self.seen_kernel_names.add(d['kernel_name'])

            return False, d
        else:
            d["reference_kernel_name"] = d["kernel_name"]
            d["reference_kernel_type"] = d["kernel_type"]


        # Generate the documentation
        if dimension == 1:
            d['edge_graph_header'] = "Edge Graph."
        elif dimension == 2:
            d['edge_graph_header'] = "X,Y Edge Graph."
        elif dimension == 3:
            d['edge_graph_header'] = "X,Y Edge Graph, sliced at z = 0."
        else:
            d['edge_graph_header'] = "X,Y Edge Graph, sliced at 0."

        ########################################
        # Now building the visualization of the graph

        x_coords = range(min(e[0] for e in edge_list) - 1,
                         max(e[0] for e in edge_list) + 2)

        y_coords = ([0] if dimension == 1 else
                    range(min(e[1] for e in edge_list) - 1,
                          max(e[1] for e in edge_list) + 2))

        node_lookup = {}

        for node in product(x_coords, y_coords):
            node_d = {}
            node_d["node_name"] = "N" + md5(str(node)).hexdigest()[:8]
            node_d["x_coord"] = graph_scale * (node[0] - min(x_coords))
            node_d["y_coord"] = graph_scale * (node[1] - min(y_coords)) if len(node) >= 2 else 0
            node_d["label"] = ("%(x_coord)d,%(y_coord)d" % node_d
                               if dimension >= 2 else node[0])
            node_d["node_color"] = "#000000"

            node_lookup[node] = node_d

        if dimension == 1:
            edge_list_2d = [y + (0,) for y in edge_list]
        elif dimension == 2:
            edge_list_2d = edge_list
        else:
            edge_list_2d = sorted(set(edge[:2] for edge in edge_list
                                      if edge[:2] != (0,0)))

        # Set all the affected node colors
        for node in edge_list_2d:
            node_lookup[node]["node_color"] = "#0000FF"

        node_lookup[(0,0)]["node_color"] = "#FF0000"

        graphviz_node_list = [T.graphviz_node % node_d
                              for k, node_d in sorted(node_lookup.iteritems())]
        graphviz_grid_edge_dict = {}
        graphviz_edge_list = []

        def edgeKey(n1, n2):
            return tuple(sorted([n1,n2]))

        # Now set the grid edges
        for node_1 in product(x_coords, y_coords):
            for xd, yd in [(0, 1), (1, 0)]:
                node_2 = (node_1[0] + xd, node_1[1] + yd)

                if node_2 in node_lookup:
                    local_edge_d = {
                        "node_1" : node_lookup[node_1]["node_name"],
                        "node_2" : node_lookup[node_2]["node_name"],
                        "line_length" : graph_scale
                    }

                    key = edgeKey(node_1, node_2)

                    graphviz_grid_edge_dict[key] = \
                        (node_1, node_2, T.graphviz_grid_edge % local_edge_d)

        # Finally, set the kernel edges
        for node in edge_list_2d:
            rk, nk = ("node_1", "node_2") if node[0] > 0 else ("node_1", "node_2")

            local_edge_d = {
                rk : node_lookup[(0,0)]["node_name"],
                nk : node_lookup[node]["node_name"],
                "line_length" : graph_scale * sqrt(node[0]**2 + node[1]**2)
            }

            key = edgeKey((0,0), node)

            graphviz_edge_list.append(T.graphviz_edge % local_edge_d)

            if key in graphviz_grid_edge_dict:
                del graphviz_grid_edge_dict[key]

        graphviz_grid_edge_list = (edge for n1,n2,edge
                                   in sorted(graphviz_grid_edge_dict.itervalues()))

        d['graphviz_nodes'] = '\n'.join(graphviz_node_list)
        d['graphviz_edges'] = '\n'.join(graphviz_edge_list)
        d['graphviz_grid_edges'] = '\n'.join(graphviz_grid_edge_list)
        d['description'] = description
        d['documentation'] = T.kernel_documentation % d
        d['sort_key'] = (dimension, radius, d["kernel_name"])

        self.seen_kernels[kernel_key] = d
        self.seen_kernel_names.add(d['kernel_name'])

        return True, d

    def addKernel(self, ktype, dimension, radius = None, description = "", template_only = False, **kwargs):
        """
        Adds the kernel to the accumulation of source files, if it is
        not already present.
        """

        new_kernel, d = self.getDictionary(ktype, dimension, radius, description, template_only, **kwargs)

        if d is None:
            return

        if not template_only:
            self.factory_kernel_dimension_list.append(
                T.factory_kernel_dimension_entry % d)

        if not new_kernel and d["reference_kernel_type"] == ktype:
            return

        ############################################################
        # Eliminate duplicate parts of the compilation

        if new_kernel:
            explicit_instantiation_filename = d['kernel_name'] + ".cpp"

            self.kernel_source_list.append(T.kernel_code % d)

            if not template_only:

                d2 = d.copy()
                explicit_instantiations = []

                for s in self.solver_list:
                    d2["solver_name"] = s

                    explicit_instantiations.append(T.explicit_instantiation_code % d2)

                    self.factory_creation_functions.append(
                        T.factory_creation_function % d2)

                    self.instantiation_suppression_list.append(
                        T.instantiation_suppression_code % d2)

                d2["explicit_instantiation_list"] = "\n".join(explicit_instantiations)

                self.instantiation_file_list.append(
                    (explicit_instantiation_filename, T.explicit_instantiation_file % d2) )

        else:
            self.kernel_source_list.append(T.kernel_duplicate_code % d)

        ############################################################
        # Common stuff

        if not template_only:

            self.kernel_list.append((d["sort_key"], d["kernel_name"]) )

            d2 = d.copy()

            for s in self.solver_list:
                d2["solver_name"] = s
                self.factory_generator_map.append(T.factory_map_entry % d2)

            self.documentation_file_list[(d['dimension'], d['documentation_file'])].append(
                (d["sort_key"], d['kernel_name'], d['documentation']))

    def write(self):
        """
        Writes everything accumulated out to the given directory.
        """

        ##################################################
        # Writing the kernel file
        self.ctx.start_msg("Writing kernels to '%s' with %d kernels." % (kernel_write_file, len(self.kernel_source_list)))

        f = open(kernel_write_file, 'w')

        kernel_d = {}
        kernel_d['kernel_source'] = '\n'.join(self.kernel_source_list)
        kernel_d['max_kernel_dimension'] = max(self.kernel_dimensions_considered)

        f.write(T.header + T.kernel_file % kernel_d)
        f.close()
        self.ctx.end_msg("done.")


        ##################################################
        # Explicit Instantiation files

        self.ctx.start_msg("Writing explicit instantiation files to %s."
                           % kernel_explicit_instantiation_directory)

        for filename, content in self.instantiation_file_list:
            f = open(join(self.kernel_explicit_instantiation_directory, filename), 'w')
            f.write(T.header + content)
            f.close()

        self.ctx.end_msg("done")

        ##################################################
        # Factory function

        self.ctx.start_msg("Writing factory function to %s." % self.factory_write_file)

        factory_d = {}
        factory_d["instantiation_suppression_list"] = (
            "\n".join(self.instantiation_suppression_list))

        factory_d["factory_creation_functions"] = (
            "\n".join(self.factory_creation_functions))

        factory_d["factory_generator_map"] = (
            ",\n".join(self.factory_generator_map))

        factory_d["factory_solver_list"] = (
            ",\n".join('"%s"' % s for s in sorted(self.solver_list)))

        factory_d["factory_kernel_list"] = (
            ",\n".join('"%s"' % s for key, s in sorted(self.kernel_list)))

        factory_d["factory_kernel_dimension_list"] = (
            ",\n".join('%s' % s for s in self.factory_kernel_dimension_list))

        f = open(self.factory_write_file, 'w')
        f.write(T.header + T.factory_file % factory_d)
        f.close()

        self.ctx.end_msg("done")


        ##################################################
        # Documentation Files

        for (dimension, filename), content_list in sorted(
            self.documentation_file_list.iteritems()):

            docfilename = join(self.documentation_directory, filename)

            self.ctx.start_msg("Writing kernel documentation to %s." % docfilename)

            f = open(docfilename, 'w')
            f.write((T.kernel_documentation_header % {"dimension" : dimension})
                    + '\n\n'.join(content for key, name, content in sorted(content_list)))
            f.close()

            self.ctx.end_msg("done")
