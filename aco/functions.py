#
# -*- coding: utf-8 -*-
#
# functions.py
#
# Auxiliary functions for smart_log_aco.py
#
# Eero Holmström, 2019
#


import parameters
import numpy as np
import copy
import random
import networkx
import smart_log_classes_aco
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.pyplot import cm
from matplotlib import colors
import matplotlib
from scipy.interpolate import CubicSpline
import itertools
import math
import sys
from PIL import Image
import warnings


#
# This is required for CSC batch jobs
#

matplotlib.use("agg")




#
# Read in a harvestability map, and slice from this map a minimal
# harvestability map containing the given bounding box. Return the
# harvestability map along with the actual bounding box of the sliced
# area.
#

def get_harvestability_map(the_graph, pile_node_sizes, pile_positions, graph_axis_limits, boundingbox_min_easting, boundingbox_max_easting, boundingbox_min_northing, boundingbox_max_northing):

        
        #
        # This is the file containing the full trafficability map
        #
        
        infile = parameters.harvestability_map_file

        
        #
        # Set this parameter so that you can read in arbitrarily large
        # .tiff files of harvestability
        #

        Image.MAX_IMAGE_PIXELS = None

        print("")
        print("Using the following parameter values:")
        print("")
        print("Upper left corner pixel center x-coordinate: %f" % parameters.harvestability_map_ulc_pixel_center_x_coordinate)
        print("Upper left corner pixel center y-coordinate: %f" % parameters.harvestability_map_ulc_pixel_center_y_coordinate)
        print("Pixel size (in m): %f" % parameters.harvestability_map_pixel_size)
        print("")
        print("Now reading in harvestability data as the image file %s..." % infile)


        #
        # Load the image
        #

        image = Image.open(infile)

        print("Done. Here are some stats on the image:")
        print("")


        #
        # Print out some stats on the image
        #

        print("Format:", image.format)
        print("Size (width, height):", image.size)
        print("Mode:", image.mode)
        print("Bands:", image.getbands())
        print("Colors (count, pixel value):", image.getcolors())
        #print("Palette:", image.getpalette())
        #print("Length of palette:", len(image.getpalette()))

        
        #
        # Create a colormap for the plot below matching the original
        # palette of the harvestability data
        #

        the_palette = image.getpalette()
        
        colors_for_pixel_values = {}

        #
        # Use black for -1
        #

        colors_for_pixel_values[-1] = (0, 0, 0)

        #
        # Get the other colors from the original palette
        #
        
        for pixel_value in [1, 2, 3, 4, 5, 6]:

                colors_for_pixel_values[pixel_value] = list(np.array(the_palette[3*pixel_value : 3*pixel_value+3]) / 255.0)

        harvestability_cmap = colors.ListedColormap([colors_for_pixel_values[key] for key in [-1, 1, 2, 3, 4, 5, 6]])
        cmap_bounds = [-1.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
        harvestability_cmap_norm = colors.BoundaryNorm(cmap_bounds, len(colors_for_pixel_values))

        
        #
        # Convert the image into a numpy array
        #

        harvestability_map = np.array(image)
        harvestability_map_y_shape = harvestability_map.shape[0]
        harvestability_map_x_shape = harvestability_map.shape[1]

        print("")
        print("Converted the image into a numpy array. Here are some stats on the array:")
        print("")
        print("Shape (y, x): (%d, %d)" % (harvestability_map_y_shape, harvestability_map_x_shape))
        print("Minimum value:", np.min(harvestability_map))
        print("Maximum value:", np.max(harvestability_map))
        print("Mean value:", np.mean(harvestability_map))


        #
        # Get the easting-northing limits of the given data
        #

        data_min_easting = parameters.harvestability_map_ulc_pixel_center_x_coordinate
        data_max_easting = parameters.harvestability_map_ulc_pixel_center_x_coordinate + parameters.harvestability_map_pixel_size*(harvestability_map_x_shape - 1)
        data_max_northing = parameters.harvestability_map_ulc_pixel_center_y_coordinate
        data_min_northing = parameters.harvestability_map_ulc_pixel_center_y_coordinate - parameters.harvestability_map_pixel_size*(harvestability_map_y_shape - 1)

        print("")
        print("The full harvestability data runs from easting %f m to %f m and northing %f m to %f m" % (data_min_easting, data_max_easting, data_min_northing, data_max_northing))
        
        
        #
        # Make sure the given bounding box is within range
        #
        
        if boundingbox_max_easting > data_max_easting or boundingbox_min_easting < data_min_easting or boundingbox_max_northing > data_max_northing or boundingbox_min_northing < data_min_northing:
                
                print("ERROR! Bounding box reaches outside of harvestability data range. Exiting.")
                exit(1)


        #
        # Slice out the harvestability data within a bounding box of the desired (harvest) area
        #
        
        #
        # Find the indeces that correspond to the closest points to
        # the given min-max boundingbox values, in such a way that the
        # bounding box fits within these indeces
        #

        i_x_max_boundingbox = (np.floor((boundingbox_max_easting - data_min_easting) / parameters.harvestability_map_pixel_size) + 1).astype(int)
        i_x_min_boundingbox = (np.floor((boundingbox_min_easting - data_min_easting) / parameters.harvestability_map_pixel_size)).astype(int)
        i_y_max_boundingbox = (np.floor(-1.0*(boundingbox_min_northing - data_max_northing) / parameters.harvestability_map_pixel_size) + 1).astype(int)
        i_y_min_boundingbox = (np.floor(-1.0*(boundingbox_max_northing - data_max_northing) / parameters.harvestability_map_pixel_size)).astype(int)

        #
        # Compute the corresponding easting and northing ranges for the harvestability map values
        #

        final_boundingbox_max_easting = parameters.harvestability_map_ulc_pixel_center_x_coordinate + i_x_max_boundingbox*parameters.harvestability_map_pixel_size
        final_boundingbox_min_easting = parameters.harvestability_map_ulc_pixel_center_x_coordinate + i_x_min_boundingbox*parameters.harvestability_map_pixel_size
        final_boundingbox_min_northing = parameters.harvestability_map_ulc_pixel_center_y_coordinate - i_y_max_boundingbox*parameters.harvestability_map_pixel_size
        final_boundingbox_max_northing = parameters.harvestability_map_ulc_pixel_center_y_coordinate - i_y_min_boundingbox*parameters.harvestability_map_pixel_size
        
        print("The given bounding box for the area runs from easting %f m to %f m and northing %f m to %f m" % (boundingbox_min_easting, boundingbox_max_easting, boundingbox_min_northing, boundingbox_max_northing))
        print("The slicing indeces for the bounding box run from i_x = %d to %d, i_y = %d to %d" % (i_x_min_boundingbox, i_x_max_boundingbox, i_y_max_boundingbox, i_y_min_boundingbox))
        print("This corresponds to a final boundingbox running from easting %f m to %f m and northing %f m to %f m" % (final_boundingbox_min_easting, final_boundingbox_max_easting, final_boundingbox_min_northing, final_boundingbox_max_northing))
        print("")

        #
        # Slice the map
        #

        harvestability_map = harvestability_map[i_y_min_boundingbox: i_y_max_boundingbox + 1, i_x_min_boundingbox: i_x_max_boundingbox + 1]
        harvestability_map_y_shape = harvestability_map.shape[0]
        harvestability_map_x_shape = harvestability_map.shape[1]

        print("Sliced the harvestability map to the given bounding box. Here are some stats:")
        print("")
        print("Shape (y, x): (%d, %d)" % (harvestability_map_y_shape, harvestability_map_x_shape))
        print("Minimum value:", np.min(harvestability_map))
        print("Maximum value:", np.max(harvestability_map))
        print("Mean value:", np.mean(harvestability_map))


        #
        # Crop out just the part of the array with values in the range
        # 1...6. Insert values of -1 elsewhere to denote missing data, which
        # can apparently be due to, e.g., bodies of water.
        #

        print("")
        print("Now setting all data cells outside of the range 1...6 to -1...")
        
        truth_mask = np.logical_and(harvestability_map >= 1, harvestability_map <= 6)
        harvestability_map = harvestability_map*(truth_mask.astype(int))
        harvestability_map[np.where(harvestability_map < 1)] = -1

        print("Done. Here are some stats on the new array:")
        print("")
        print("Minimum value:", np.min(harvestability_map))
        print("Maximum value:", np.max(harvestability_map))
        print("Mean value:", np.mean(harvestability_map))
        print("Standard deviation:", np.std(harvestability_map))
        print("")

        hist, bin_edges = np.histogram(harvestability_map, range = [1, 7], bins = 6)

        print("Histogram of values for harvestability of 1, 2, 3, 4, 5, 6:")
        print("")
        print(hist)
        #print(bin_edges)
        print("")

        #
        # Plot the harvestability map
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Harvestability')
        
        #
        # Draw the full graph into the background first, along with
        # the road connectivity
        #                

        networkx.draw(the_graph, node_color = 'white', node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)


        #
        # Then, draw the harvestability map. Use the same color
        # palette as Metsäkeskus uses.
        #
        
        plt.imshow(harvestability_map, extent = (final_boundingbox_min_easting - parameters.harvestability_map_pixel_size/2.0, final_boundingbox_max_easting + parameters.harvestability_map_pixel_size/2.0, final_boundingbox_min_northing - parameters.harvestability_map_pixel_size/2.0, final_boundingbox_max_northing + parameters.harvestability_map_pixel_size/2.0), cmap = harvestability_cmap, norm = harvestability_cmap_norm)
        
        plt.axis(graph_axis_limits)
        plt.gca().set_aspect('equal', adjustable = 'box')
        cbar = plt.colorbar(ticks = [1, 2, 3, 4, 5, 6], pad = 0.06, shrink = 0.55)

        if parameters.save_figures_to_files:
                        
                plt.savefig('harvestability.png')

                
        plt.close()
        
        return harvestability_map, final_boundingbox_min_easting, final_boundingbox_max_easting, final_boundingbox_min_northing, final_boundingbox_max_northing


#
# Decompose the full, real-life graph, into smaller "region graphs" for which the
# optimization will be performed. Do this through the following steps:
#
# 1. Find the road strips between intersections in the real-life
# graph.
#
# 2. Get fragments remaining outside of these road strips, i.e., the
# intersection nodes.
#
# 3. Join each fragment to a road strip. The decomposition will be
# built from the resulting set of road strips, some of which have been
# augmented with fragments. We collectively call these "subgraphs".
#
# 4. Enforce a minimum size on the subgraphs by merging subgraphs
# together until all subgraphs are of at least the minimum size.
#
# 5. Build the decomposition by first choosing the closest subgraph to
# the roadside node that is not already in the decomposition, then
# adding subgraphs in order of increasing distance from this "seed
# subgraph", keep on adding subgraphs until the region graph is at
# least of the minimum required size.
#
# Return the decomposition, which is a list of lists of node indeces,
# each list giving the nodes for that region graph.
#

def decompose_graph(the_original_graph, original_pile_node_colors, original_pile_node_sizes, pile_color_legend_handles, positions, graph_axis_limits):

        #
        # Ignore (deprecation) warnings for this function
        #
        
        warnings.simplefilter('ignore') 

        #
        # Remove the roadside node from the graph and from the node
        # size and node color lists. We deal here only with the pickup
        # nodes.
        #
        
        the_graph = copy.deepcopy(the_original_graph)
        pile_node_colors = copy.deepcopy(original_pile_node_colors)
        pile_node_sizes = copy.deepcopy(original_pile_node_sizes)
        
        the_graph.remove_node(0)
        pile_node_colors.pop(0)
        pile_node_sizes.pop(0)

        
        #
        # This will be useful later on
        #
        
        list_of_all_nodes = list(the_graph.nodes)

        
        #
        # 1. Find the road strips in the real-life graph
        #

        print("")
        print("Now finding the road strips in the real-life graph...")

        
        #
        # First, find all nodes of degree >= 3, i.e., the intersections
        #

        intersection_nodes = [node[0] for node in the_graph.degree(list(the_graph.nodes)) if node[1] >= 3]

        
        #
        # Then, remove these intersection nodes from the graph to find the road strips
        #

        graph_of_all_road_strips = copy.deepcopy(the_graph)
        
        graph_of_all_road_strips.remove_nodes_from(intersection_nodes)

        list_of_road_strip_graphs = [graph_of_all_road_strips.subgraph(this_road_strip_nodes).copy() for this_road_strip_nodes in networkx.connected_components(graph_of_all_road_strips)]
        
        road_strip_node_lists = [list(this_road_strip_nodes) for this_road_strip_nodes in networkx.connected_components(graph_of_all_road_strips)]
        
        number_of_road_strips = len(road_strip_node_lists)
        

        #
        # This will be useful later on
        #
        
        nodes_in_road_strips = []
        
        for this_road_strip in road_strip_node_lists:
                
                nodes_in_road_strips.extend(this_road_strip)

        nodes_in_road_strips = set(nodes_in_road_strips)
        
        
        #
        # Compute and print out some statistics on the road strips
        #

        road_strip_lengths = []
        
        for this_road_strip in road_strip_node_lists:

                road_strip_lengths.append(len(this_road_strip))
        
        road_strip_lengths = np.array(road_strip_lengths)

        minimum_road_strip_length = np.min(road_strip_lengths)
        mean_road_strip_length = np.mean(road_strip_lengths)
        maximum_road_strip_length = np.max(road_strip_lengths)
        std_of_road_strip_length = np.std(road_strip_lengths)
        
        print("Done. Extracted a total of %d road strips. Road strip size ranges from %d to %d nodes, with a mean of %f and a standard deviation of %f" % (number_of_road_strips, minimum_road_strip_length, maximum_road_strip_length, mean_road_strip_length, std_of_road_strip_length))

        
        #
        # Plot the road strips on top of the full real-life graph, using a
        # different color for each road strip
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Road strips')

        
        #
        # Draw the full graph into the background first, along with
        # the road connectivity
        #                

        networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)


        #
        # Then, draw all nodes belonging to road strips in gray
        #
                
        the_road_strip_graph = copy.deepcopy(the_graph)
        all_nodes_except_those_in_road_strips = set(list_of_all_nodes).difference(nodes_in_road_strips)
        the_road_strip_graph.remove_nodes_from(all_nodes_except_those_in_road_strips)
                                                       
        networkx.draw(the_road_strip_graph, node_color = 'grey', alpha = 0.5, node_size = parameters.node_size_in_road_strip_graph, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 's', pos = positions, edgelist = [], with_labels = parameters.draw_graph_nodes_with_labels)

        
        #
        # Create a color space for the edges of the road strips
        #
        
        road_strip_colors = cm.rainbow(np.linspace(0, 1, number_of_road_strips))
        
        #
        # Shuffle the list of road strip colors, so that the different
        # road strips are easier to make out in the plot
        #

        random.shuffle(road_strip_colors)
        
        i_road_strip = 0

        for this_road_strip_graph in list_of_road_strip_graphs:
                
                #
                # Get the path edges for plotting this road strip
                #
                
                this_road_strip_edges = this_road_strip_graph.edges

                
                #
                # Create a list of colors
                #

                this_road_strip_edge_colors = [road_strip_colors[i_road_strip] for i in range(0, len(this_road_strip_edges))]

                
                #
                # Then draw the path for this road strip
                #
                
                networkx.draw_networkx_edges(this_road_strip_graph, pos = positions, edgelist = this_road_strip_edges, edge_color = this_road_strip_edge_colors, width = parameters.road_strip_width)
                
                i_road_strip = i_road_strip + 1

                
        plt.legend(handles = pile_color_legend_handles, loc = 'upper right')
        plt.axis(graph_axis_limits)
        plt.gca().set_aspect('equal', adjustable = 'box')
    
        if parameters.save_figures_to_files:
                        
                plt.savefig('road_strips_of_real-life_graph.png')



        plt.close()
        
        
        #
        # 2. Find all nodes outside of the road strips, i.e., the fragments
        # left over from the road strips, i.e., the intersection nodes
        #

        print("")
        print("Now finding the fragments outside of the road strips...")

        
        #
        # Create a new graph object for this purpose
        #

        graph_of_all_fragments = copy.deepcopy(the_graph)

        
        #
        # Remove the nodes belonging to road strips from the graph
        #
        
        graph_of_all_fragments.remove_nodes_from(nodes_in_road_strips)


        #
        # Get each fragment into a separate graph. Get the fragments
        # also into a simple list of node lists.
        #

        list_of_fragment_graphs = [graph_of_all_fragments.subgraph(this_component_nodes).copy() for this_component_nodes in networkx.connected_components(graph_of_all_fragments)]

        fragment_node_lists = [list(this_fragment_graph) for this_fragment_graph in networkx.connected_components(graph_of_all_fragments)]

        
        #
        # Print out some statistics on the fragments
        #
        
        number_of_fragments = len(list_of_fragment_graphs)

        fragment_lengths = []
        
        if number_of_fragments > 0:
        
                for this_fragment in fragment_node_lists:

                        fragment_lengths.append(len(this_fragment))
                        
        else:

                fragment_lengths = [0]

                
        fragment_lengths = np.array(fragment_lengths)
        
        minimum_size_of_fragment = np.min(fragment_lengths)
        maximum_size_of_fragment = np.max(fragment_lengths)
        mean_size_of_fragment = np.mean(fragment_lengths)
        std_of_fragment_size = np.std(fragment_lengths)

        print("Done. Found a total of %d fragments. Fragment size ranges from %d to %d nodes, with a mean of %f and a standard deviation of %f" % (number_of_fragments, minimum_size_of_fragment, maximum_size_of_fragment, mean_size_of_fragment, std_of_fragment_size))

                
        #
        # Plot the fragments
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Fragments')

        #
        # Draw the full graph into the background first, along with
        # the road connectivity
        #                

        networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)        

        #
        # Then, plot the nodes in the fragments
        #
                
        networkx.draw(graph_of_all_fragments, node_color = 'grey', node_size = parameters.fragment_graph_node_size, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'd', pos = positions, edgelist = [], with_labels = parameters.draw_graph_nodes_with_labels)


        #
        # Then, plot the edges of each fragment in a distinct color
        #
        
        i_fragment = 0
        
        fragment_colors = cm.rainbow(np.linspace(0, 1, number_of_fragments))

        for this_fragment_graph in list_of_fragment_graphs:

                #
                # Get the path edges for plotting this fragment
                #

                this_fragment_edges = this_fragment_graph.edges

                
                #
                # Create a list of colors
                #

                this_fragment_edge_colors = [fragment_colors[i_fragment] for i in range(0, len(this_fragment_edges))]

                
                #
                # Then draw the path for this fragment
                #
                
                networkx.draw_networkx_edges(this_fragment_graph, pos = positions, edgelist = this_fragment_edges, edge_color = this_fragment_edge_colors, width = parameters.fragment_path_width)
                
                i_fragment = i_fragment + 1


        plt.legend(handles = pile_color_legend_handles, loc = 'upper right')
        plt.axis(graph_axis_limits)
        plt.gca().set_aspect('equal', adjustable = 'box')
    
        if parameters.save_figures_to_files:
                        
                plt.savefig('fragments_for_real-life_graph.png')

        plt.close()

        
        #
        # 3. Join each fragment to the first road strip that it is
        # directly connected to in the real-life graph. The set of
        # subgraphs used to build the decomposition will be made up of
        # these augmented road strips as well as any remaining road
        # strips.
        #
        
        print("")
        print("Now joining each fragment to the first road strip that it is directly connected to...")

        
        list_of_subgraphs = []
        list_of_road_strips_in_subgraphs = []

        n_joining_operations = 0
        
        for this_fragment_node_list in fragment_node_lists:

                #
                # Simply see if this fragment and this road strip have any two nodes,
                # respectively, which are connected by an edge in the real-life graph
                #

                #
                # Get a list of all nodes that are involved in edges involved with
                # the nodes of this fragment
                #
                
                edges_of_this_fragment_nodes = the_graph.edges(this_fragment_node_list)
                nodes_in_edges_of_this_fragment_0 = [e[0] for e in edges_of_this_fragment_nodes]
                nodes_in_edges_of_this_fragment_1 = [e[1] for e in edges_of_this_fragment_nodes]
                
                nodes_in_edges_of_this_fragment = []
                nodes_in_edges_of_this_fragment.extend(nodes_in_edges_of_this_fragment_0)
                nodes_in_edges_of_this_fragment.extend(nodes_in_edges_of_this_fragment_1)
                
                #
                # BEGIN DEBUG
                #

                #print("")
                #print("fragment node list:", this_fragment_node_list)
                #print("the edge view of this fragment:", the_graph.edges(this_fragment_node_list))
                #print("edge nodes of this fragment:", nodes_in_edges_of_this_fragment)
                
                #
                # END DEBUG
                #

                i_road_strip = 0
                
                for this_road_strip in road_strip_node_lists:


                        #
                        # BEGIN DEBUG
                        #
                        
                        #print("nodes of this road strip:", this_road_strip)
                        
                        #
                        # END DEBUG
                        #
                        
                        
                        #
                        # See if this fragment is connected to this road strip
                        #

                        if len(set(nodes_in_edges_of_this_fragment).intersection(set(this_road_strip))) > 0:

                                #
                                # BEGIN DEBUG
                                #

                                #print("SUCCESS! Now joining fragment with this road strip.")
                                
                                #
                                # END DEBUG
                                #
                                
                                #
                                # You have found a road strip that is connected to this fragment. Create
                                # the joint subgraph.
                                #
                                
                                this_joint_subgraph = copy.deepcopy(the_graph)

                                nodes_in_this_road_strip = set(this_road_strip)
                                nodes_in_this_fragment = set(this_fragment_node_list)
                                nodes_in_this_subgraph = nodes_in_this_road_strip.union(nodes_in_this_fragment)
                
                                all_nodes_except_those_in_this_subgraph = set(list_of_all_nodes).difference(nodes_in_this_subgraph)
                        
                                this_joint_subgraph.remove_nodes_from(all_nodes_except_those_in_this_subgraph)

                                list_of_subgraphs.append(this_joint_subgraph)

                                n_joining_operations = n_joining_operations + 1

                                
                                #
                                # Also, keep track of which road strips are included in the subgraphs already
                                #
                                                        
                                list_of_road_strips_in_subgraphs.append(i_road_strip)

                                
                                #
                                # Plot the road strip and fragment before and after joining
                                #

                                #
                                #  Before
                                #

                                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                                plt.rc('font', size = parameters.default_font_size)
                                plt.title('Road strip (red) and fragment (blue)')

                                #
                                # Draw the full graph into the background first, along with
                                # the road connectivity
                                #                
        
                                #networkx.draw(the_graph, node_color = 'black', node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                
                                this_road_strip_graph = copy.deepcopy(the_graph)
                                all_nodes_except_those_in_road_strip = set(list_of_all_nodes).difference(nodes_in_this_road_strip)
                                this_road_strip_graph.remove_nodes_from(all_nodes_except_those_in_road_strip)
                                                       
                                networkx.draw(this_road_strip_graph, node_color = 'red', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                        
                                this_fragment_graph = copy.deepcopy(the_graph)
                                all_nodes_except_those_in_fragment = set(list_of_all_nodes).difference(nodes_in_this_fragment)
                                this_fragment_graph.remove_nodes_from(all_nodes_except_those_in_fragment)
                                                        
                                networkx.draw(this_fragment_graph, node_color = 'blue', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

                                plt.axis(graph_axis_limits)

                                if parameters.save_figures_to_files and parameters.save_joining_and_merging_figures_to_files:

                                        plt.savefig('joining_operation_' + str(n_joining_operations) + '_before.png')

                                plt.close()
                                

                                #
                                # After
                                #
                                                        
                                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                                plt.rc('font', size = parameters.default_font_size)
                                plt.title('Joined')

                                #
                                # Draw the full graph into the background first, along with
                                # the road connectivity
                                #                
        
                                #networkx.draw(the_graph, node_color = 'black', node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                
                                networkx.draw(this_joint_subgraph, node_color = 'green', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                
                                plt.axis(graph_axis_limits)
                                
                                if parameters.save_figures_to_files and parameters.save_joining_and_merging_figures_to_files:

                                        plt.savefig('joining_operation_' + str(n_joining_operations) + '_after.png')

                                plt.close()
                                
                                break

                        i_road_strip = i_road_strip + 1        
        

        print("Done. Performed a total of %d joining operations." % n_joining_operations)

        
        #
        # Add the road strips which are not yet parts of subgraphs into the
        # list of subgraphs.
        #
        
        i_road_strip = 0
        
        for this_road_strip in road_strip_node_lists:

                #
                # DEBUG
                #
                
                #print("This road strip:", this_road_strip)

                #
                # END DEBUG
                #
                
                if i_road_strip not in list_of_road_strips_in_subgraphs:

                        #
                        # Form a graph object of this road strip and append it to the list of subgraphs
                        #

                        graph_of_this_road_strip = copy.deepcopy(the_graph)
                                
                        nodes_in_this_road_strip = set(this_road_strip)
                        
                        all_nodes_except_those_in_this_road_strip = set(list_of_all_nodes).difference(nodes_in_this_road_strip)
                                
                        graph_of_this_road_strip.remove_nodes_from(all_nodes_except_those_in_this_road_strip)
                                                                
                        list_of_subgraphs.append(graph_of_this_road_strip)

                        
                i_road_strip = i_road_strip + 1




        #
        # 4. Enforce a minimum size on the subgraphs. If subgraph A is
        # smaller than the minimum required size, merge subgraph A
        # with another subgraph B. The protocol for finding subgraph B
        # is the following:
        #
        # I) Attempt to find a subgraph B which has nodes in common
        # with subgraph A. If this fails, go to II).
        #
        # II) Attempt to find a subgraph B which is directly connected
        # to subgraph A in the real-life graph.
        #
        # Loop over the subgraphs until all of them are at least as large
        # as the minimum required size.
        #

        print("")
        print("Now checking that all subgraphs are larger than the minimum required size, and merging subgraphs if necessary...")

        
        all_subgraphs_larger_than_minimum_size = False

        n_merging_operations = 0
        
        
        while not all_subgraphs_larger_than_minimum_size:

                #
                # BEGIN DEBUG
                #

                #print("Starting new loop.")
                
                #
                # END DEBUG
                #
                
                all_subgraphs_larger_than_minimum_size = True
                
                for this_subgraph in list_of_subgraphs:
                        
                        if len(this_subgraph.nodes) < parameters.minimum_number_of_nodes_per_subgraph:

                                all_subgraphs_larger_than_minimum_size = False

                                subgraph_to_merge_with = None
                        
                                #
                                # I) Try to find a subgraph that shares one or more node with this subgraph
                                #

                                for other_subgraph in list_of_subgraphs:

                                        if this_subgraph is not other_subgraph:

                                                #
                                                # BEGIN DEBUG
                                                #
                                                
                                                #print("The two graphs, see if they share nodes:", set(list(this_subgraph.nodes)), set(list(other_subgraph.nodes)))
                                                
                                                #
                                                # END DEBUG
                                                #

                                                
                                                if len(set(list(this_subgraph.nodes)).intersection(set(list(other_subgraph.nodes)))) > 0:

                                                        #
                                                        # Success! Merge the two subgraphs.
                                                        #
                                                        
                                                        subgraph_to_merge_with = other_subgraph

                                                        merged_subgraph = networkx.compose(this_subgraph, subgraph_to_merge_with)

                                                        n_merging_operations = n_merging_operations + 1
                                                        
                                                        #
                                                        # DEBUG
                                                        #

                                                        #print("SUCCESS (1)! Found subgraph B with common nodes with A.")
                                                        
                                                        #
                                                        # END DEBUG
                                                        #
                                                        
                                                        break

                                

                                #
                                # II) If we were not successful, attempt to find a subgraph which is
                                # directly connected to this subgraph
                                #

                                if subgraph_to_merge_with is None:

                                        edges_of_this_subgraph_nodes = the_graph.edges(list(this_subgraph.nodes))
                                        nodes_in_edges_of_this_subgraph_0 = [e[0] for e in edges_of_this_subgraph_nodes]
                                        nodes_in_edges_of_this_subgraph_1 = [e[1] for e in edges_of_this_subgraph_nodes]
                
                                        nodes_in_edges_of_this_subgraph = []
                                        nodes_in_edges_of_this_subgraph.extend(nodes_in_edges_of_this_subgraph_0)
                                        nodes_in_edges_of_this_subgraph.extend(nodes_in_edges_of_this_subgraph_1)
                                        
                                        for other_subgraph in list_of_subgraphs:

                                                if this_subgraph is not other_subgraph:

                                                        #
                                                        # BEGIN DEBUG
                                                        #
                                                        
                                                        #print("The two graphs, see if they are connected:", set(list(this_subgraph.nodes)), set(list(other_subgraph.nodes)))

                                                        #
                                                        # END DEBUG
                                                        #
                                                        
                                                        nodes_of_other_subgraph = list(other_subgraph.nodes)

                                                        if len(set(nodes_in_edges_of_this_subgraph).intersection(set(nodes_of_other_subgraph))) > 0:
                                                        
                                                                #
                                                                # Success! Merge the two subgraphs.
                                                                #
                                                        
                                                                subgraph_to_merge_with = other_subgraph

                                                                merged_subgraph = copy.deepcopy(the_graph)
                                                                
                                                                nodes_in_this_subgraph = set(list(this_subgraph.nodes))
                                                                nodes_in_other_subgraph = set(list(other_subgraph.nodes))
                                                                nodes_in_merged_subgraph = nodes_in_this_subgraph.union(nodes_in_other_subgraph)
                                                                
                                                                all_nodes_except_those_in_this_subgraph = set(list_of_all_nodes).difference(nodes_in_merged_subgraph)
                                                                
                                                                merged_subgraph.remove_nodes_from(all_nodes_except_those_in_this_subgraph)
                                                                                                                                
                                                                n_merging_operations = n_merging_operations + 1
                                                                
                                                                #
                                                                # DEBUG
                                                                #

                                                                #print("SUCCESS (2)! Found subgraph B directly connected to A.")
                                                                
                                                                #
                                                                # END DEBUG
                                                                #

                                                                break
                                
                                        
                                if subgraph_to_merge_with is None:

                                        print("ERROR! Something went wrong with merging subgraphs. Exiting.")
                                        exit(1)

                                        
                                #
                                # Remove the two separate subgraphs from the list of subgraphs
                                #

                                i_sg = 0
                                
                                for sg in list_of_subgraphs:

                                        if sg is this_subgraph:

                                                i_remove = i_sg
                                                break
                                                        
                                        i_sg = i_sg + 1
                                
                                list_of_subgraphs.pop(i_remove)

                                i_sg = 0
                                                        
                                for sg in list_of_subgraphs:

                                        if sg is subgraph_to_merge_with:

                                                i_remove = i_sg
                                                break

                                        i_sg = i_sg + 1
                                
                                list_of_subgraphs.pop(i_remove)
                                
                                
                                #
                                # Add the newly formed subgraph to the list of subgraphs
                                #

                                list_of_subgraphs.append(merged_subgraph)


                                #
                                # Plot the two subgraphs before and after merging
                                #

                                #
                                # Before
                                #
                                
                                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                                plt.rc('font', size = parameters.default_font_size)
                                plt.title('Subgraph A (red) and subgraph B (blue)')

                                
                                #
                                # Draw the full graph into the background first, along with
                                # the road connectivity
                                #                
        
                                #networkx.draw(the_graph, node_color = 'black', node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                
                                networkx.draw(this_subgraph, node_color = 'red', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                                        
                                networkx.draw(subgraph_to_merge_with, node_color = 'blue', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

                                plt.axis(graph_axis_limits)
                                
                                if parameters.save_figures_to_files and parameters.save_joining_and_merging_figures_to_files:

                                        plt.savefig('merging_operation_' + str(n_merging_operations) + '_before.png')

                                plt.close()


                                #
                                # After
                                #
                                                        
                                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                                plt.rc('font', size = parameters.default_font_size)
                                plt.title('Merged')

                                #
                                # Draw the full graph into the background first, along with
                                # the road connectivity
                                #                
        
                                #networkx.draw(the_graph, node_color = 'black', node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)                                
                                
                                networkx.draw(merged_subgraph, node_color = 'green', alpha = 0.5, node_size = parameters.node_size_in_before_and_after_joining_plots, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
                                
                                plt.axis(graph_axis_limits)

                                if parameters.save_figures_to_files and parameters.save_joining_and_merging_figures_to_files:

                                        plt.savefig('merging_operation_' + str(n_merging_operations) + '_after.png')

                                plt.close()

                                break


        print("Done. Performed a total of %d merging operations." % n_merging_operations)

        
        #
        # Now the full list of subgraphs has been constructed. Check
        # that each node in each subgraph can be accessed from any
        # other node. This is to make sure you get nice and smooth
        # coverage when creating the decomposition later on.
        #

        number_of_subgraphs = len(list_of_subgraphs)
        
        
        print("")
        print("Now checking that in each subgraph, each node can be accessed from any other node...")
        
        
        i_subgraph = 0
        
        for this_subgraph in list_of_subgraphs:
                
                if not networkx.is_connected(this_subgraph):
                                                                                
                        print("ERROR! Subgraph of index %d is not connected. Exiting." % i_subgraph)
                        print("Here's the subgraph:", list(this_subgraph.nodes), this_subgraph.edges)
                        exit(1)
        
                i_subgraph = i_subgraph + 1
        
        
        print("Done.")
        
                
        #
        # Also, check that each node is included at least once in the
        # ensemble of subgraphs.
        #

        #
        # Find set of nodes included in the subgraphs
        #
        
        nodes_in_subgraphs = []
        
        for this_subgraph in list_of_subgraphs:

                nodes_in_this_subgraph = list(this_subgraph.nodes)
                nodes_in_subgraphs.extend(nodes_in_this_subgraph)
        
        nodes_in_subgraphs = set(nodes_in_subgraphs)

        
        #
        # Then, see if any nodes are missing
        #
        
        for node in list_of_all_nodes:

                if node not in nodes_in_subgraphs:

                        print("ERROR! Node %d missing from subgraphs. Exiting." % node)
                        exit(1)


        #
        # Compute and print out some stats on the subgraphs
        #
        
        subgraph_lengths = []
        
        for this_subgraph in list_of_subgraphs:

                subgraph_lengths.append(len(this_subgraph.nodes))

        subgraph_lengths = np.array(subgraph_lengths)
        
        minimum_size_of_subgraph = np.min(subgraph_lengths)
        maximum_size_of_subgraph = np.max(subgraph_lengths)
        mean_size_of_subgraph = np.mean(subgraph_lengths)
        std_of_subgraph_size = np.std(subgraph_lengths)

        print("")
        print("Formed a total of %d subgraphs. Subgraph size ranges from %d to %d, with a mean of %f and a standard deviation of %f" % (number_of_subgraphs, minimum_size_of_subgraph, maximum_size_of_subgraph, mean_size_of_subgraph, std_of_subgraph_size))


        #
        # So now you have a set of subgraphs which cover the entire,
        # full, real-life graph. Each subgraph is connected.
        #
        # Next, go ahead and create the decomposition of the real-life
        # graph. Do this through the following approach:
        #
        # 1. Order the subgraphs by increasing distance of the center
        #    of mass of each from the roadside node.
        #
        # 2. Initialize a new region graph in the decomposition. Each
        #    region graph is described only by the set of nodes it
        #    comprises. Start by adding the subgraph closest to the
        #    roadside node. We call this "subgraph A". Remove subgraph
        #    A from the list of available subgraphs.
        #
        # 3. Continue adding subgraphs in order of increasing distance
        #    from subgraph A. When adding a subgraph, subtract from it
        #    the nodes already in the region graph.
        #
        # 4. Stop when you pass the threshold for the allowed maximum
        #    number of nodes in a decomposition region graph.
        #
        # 5. If nodes still remain outside of the region graphs, go
        #    back to part 2.
        #


        #
        # Compute the center of mass for each subgraph
        #

        cms_of_the_subgraphs = []
        
        for this_subgraph in list_of_subgraphs:

                this_cm = np.array([0.0, 0.0])
                this_number_of_nodes = 0
                
                for node in this_subgraph.nodes:

                        this_number_of_nodes = this_number_of_nodes + 1

                        this_node_x = this_subgraph.nodes[node]['data'].x_position
                        this_node_y = this_subgraph.nodes[node]['data'].y_position
                        
                        this_cm = this_cm + np.array([this_node_x, this_node_y])

                this_cm = this_cm / this_number_of_nodes

                cms_of_the_subgraphs.append(this_cm)        
        
        cms_of_the_subgraphs = np.array(cms_of_the_subgraphs)


        #
        # Plot the subgraphs on top of the real-life graph, first all into the same plot
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Subgraphs')
        
        #
        # Draw the full graph into the background first, along with
        # the road connectivity
        #                
        
        networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)
        
        #
        # Create a color space for the edges
        #
        
        subgraph_colors = cm.rainbow(np.linspace(0, 1, number_of_subgraphs))
        
        #
        # Shuffle the list of subgraph colors, so that the different
        # subgraphs are easier to make out in the plot
        #

        random.shuffle(subgraph_colors)

        #
        # Then, draw the subgraphs
        #
        
        i_subgraph = 0
        
        for this_subgraph in list_of_subgraphs:

                
                this_subgraph_nodes = list(this_subgraph.nodes)
               
                #
                # Get the path edges for plotting this subgraph
                #
        
                this_subgraph_edges = this_subgraph.edges

                #
                # Create a list of colors, for coloring each edge in
                # this subraph with the same color
                #
                
                this_subgraph_edge_colors = [subgraph_colors[i_subgraph] for i in range(0, len(this_subgraph_edges))]

                #
                # Then draw the edges for this subgraph
                #
                
                networkx.draw_networkx_edges(the_graph, pos = positions, edgelist = this_subgraph_edges, edge_color = this_subgraph_edge_colors, width = parameters.subgraph_path_width)

                plt.legend(handles = pile_color_legend_handles, loc = 'upper right')        
                plt.axis(graph_axis_limits)
                plt.gca().set_aspect('equal', adjustable = 'box')

                i_subgraph = i_subgraph + 1

                
        if parameters.save_figures_to_files:
        
                plt.savefig('subgraphs.png')

        plt.close()

        
        #
        # Plot the subgraphs on top of the full real-life graph, each in a separate figure
        #
        
        i_subgraph = 0
        
        for this_subgraph in list_of_subgraphs:
                
                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                plt.rc('font', size = parameters.default_font_size)
                plt.title('Subgraph of index ' + str(i_subgraph))
        
                #
                # Draw the full graph into the background first, along with
                # the road connectivity
                #                
        
                networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

                
                this_subgraph_nodes = list(this_subgraph.nodes)
               
                #
                # Get the path edges for plotting this subgraph
                #
        
                this_subgraph_edges = this_subgraph.edges

                #
                # Create a list of colors, for coloring each edge in
                # this subraph with the same color
                #
                
                this_subgraph_edge_colors = [subgraph_colors[i_subgraph] for i in range(0, len(this_subgraph_edges))]

                #
                # Then draw the edges for this subgraph
                #
                
                networkx.draw_networkx_edges(the_graph, pos = positions, edgelist = this_subgraph_edges, edge_color = this_subgraph_edge_colors, width = parameters.subgraph_path_width)

                plt.legend(handles = pile_color_legend_handles, loc = 'upper right')        
                plt.axis(graph_axis_limits)
                plt.gca().set_aspect('equal', adjustable = 'box')
        
                if parameters.save_figures_to_files and parameters.save_individual_subgraph_plots_to_files:
        
                        plt.savefig('subgraph_' + str(i_subgraph) + '.png')

                plt.close()

                i_subgraph = i_subgraph + 1
                                

                
        #
        # Plot the subgraphs with their centers of mass on top of the real-life graph
        #
        
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Subgraphs with centers of mass')
        
        #
        # Draw the full graph into the background first, along with the road connectivity
        #                
        
        networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

        #
        # Then, draw the subgraphs
        #
        
        i_subgraph = 0
        
        for this_subgraph in list_of_subgraphs:
        
                this_subgraph_nodes = list(this_subgraph.nodes)
               
                #
                # Get the path edges for plotting this subgraph
                #
        
                this_subgraph_edges = this_subgraph.edges

                #
                # Create a list of colors, for coloring each edge in
                # this subraph with the same color
                #
                
                this_subgraph_edge_colors = [subgraph_colors[i_subgraph] for i in range(0, len(this_subgraph_edges))]

                #
                # Then draw the edges for this subgraph
                #
                
                networkx.draw_networkx_edges(the_graph, pos = positions, edgelist = this_subgraph_edges, edge_color = this_subgraph_edge_colors, width = parameters.subgraph_path_width)

                #
                # Draw the center of mass for this subgraph
                #
                
                plt.plot(cms_of_the_subgraphs[i_subgraph][0], cms_of_the_subgraphs[i_subgraph][1], 'x', color = subgraph_colors[i_subgraph], markersize = parameters.subgraph_cm_marker_size, markeredgewidth = parameters.subgraph_cm_marker_edge_width)

                #
                # Draw the subgraph index also
                #

                plt.text(cms_of_the_subgraphs[i_subgraph][0], cms_of_the_subgraphs[i_subgraph][1], i_subgraph)
                
                plt.legend(handles = pile_color_legend_handles, loc = 'upper right')        
                plt.axis(graph_axis_limits)
                plt.gca().set_aspect('equal', adjustable = 'box')              

                i_subgraph = i_subgraph + 1

                
        if parameters.save_figures_to_files:
        
                plt.savefig('subgraphs_with_centers_of_mass.png')

        plt.close()
        
        
        #
        # Create the decomposition
        #

        print("")
        print("Now creating the decomposition from these subgraphs...")
        

        #
        # This list holds the final decomposition of the real-life
        # graph. That is, the lists of nodes each for which to perform
        # the route optimization separately.
        #
        
        decomposition = []

        #
        # These are auxiliary objects
        #
        
        nodes_in_decomposition = []
                
        nodes_remaining = copy.deepcopy(list_of_all_nodes)
        
        
        #
        # Keep creating region graphs until all nodes have been
        # included in the decomposition
        #


        #
        # Sort the subgraphs by increasing distance from the roadside node
        #
        
        roadside_node_position = np.array([parameters.x_roadside_position, parameters.y_roadside_position])

        subgraph_distances_from_roadside_node = np.linalg.norm(cms_of_the_subgraphs - roadside_node_position, axis = 1)

        subgraph_indeces_sorted_by_distance_from_roadside_node = np.argsort(subgraph_distances_from_roadside_node)

        #
        # BEGIN DEBUG
        #

        #print("")
        #print("Start: subgraph distances from roadside node:", subgraph_distances_from_roadside_node)
        #print("Start: the sorted indeces:", subgraph_indeces_sorted_by_distance_from_roadside_node)
        #print("")
        
        #
        # END DEBUG
        #

        #
        # Keep track of which subgraphs you have
        # already used up for this decomposition
        #
                
        indeces_of_subgraphs_used_for_this_decomposition = []

        
        while len(nodes_remaining) > 0:

                
                this_region_graph = []

                
                #
                # First, add the closest subgraph to the roadside node that
                # is still available. We call this the "seed subgraph".
                #

                i_subgraphs_sorted_by_distance = 0
                
                i_seed_subgraph = subgraph_indeces_sorted_by_distance_from_roadside_node[i_subgraphs_sorted_by_distance]

                #
                # BEGIN DEBUG
                #

                #print("")
                #print("indeces_of_subgraphs_used_for_this_decomposition:", indeces_of_subgraphs_used_for_this_decomposition)
                #print("")
                
                #
                # END DEBUG
                #
                
                while i_seed_subgraph in indeces_of_subgraphs_used_for_this_decomposition:

                        i_subgraphs_sorted_by_distance = i_subgraphs_sorted_by_distance + 1
                        i_seed_subgraph = subgraph_indeces_sorted_by_distance_from_roadside_node[i_subgraphs_sorted_by_distance]

                
                the_seed_subgraph = list_of_subgraphs[i_seed_subgraph]

                the_seed_subgraph_nodes = list(the_seed_subgraph.nodes)

                #
                # BEGIN DEBUG
                #

                #print("")
                #print("Seed subgraph index:", i_seed_subgraph)
                #print("")
                
                #
                # END DEBUG
                #
                
                
                #
                # Remember to first subtract out the nodes
                # that are already in the decomposition
                #

                the_seed_subgraph_nodes_to_add = list(set(the_seed_subgraph_nodes).difference(set(nodes_in_decomposition)))                        

                
                #
                # Then, add the subgraph to the region graph
                #
                        
                this_region_graph.extend(the_seed_subgraph_nodes_to_add)
                indeces_of_subgraphs_used_for_this_decomposition.append(i_seed_subgraph)

                
                #
                # Keep track of which nodes are already in the
                # decomposition
                #
                        
                nodes_in_decomposition.extend(the_seed_subgraph_nodes_to_add)

                
                #
                # Remove the nodes that you've used for this
                # region graph from the list of nodes that still
                # need to be added
                #

                nodes_remaining = list(set(nodes_remaining).difference(set(the_seed_subgraph_nodes_to_add)))


                #
                # Then, order the remaining subgraphs in terms of
                # distance from the center of mass of the seed
                # subgraph
                #

                seed_subgraph_cm = cms_of_the_subgraphs[i_seed_subgraph, :]

                distances_of_other_subgraphs_from_the_seed_subgraph = np.linalg.norm(cms_of_the_subgraphs - seed_subgraph_cm, axis = 1)

                subgraph_indeces_sorted_by_distance_from_the_seed_subgraph = np.argsort(distances_of_other_subgraphs_from_the_seed_subgraph)
                

                #
                # BEGIN DEBUG
                #

                #print("")
                #print("the sorted indeces by distance from seed subgraph:", subgraph_indeces_sorted_by_distance_from_the_seed_subgraph)
                #print("")
                
                #
                # END DEBUG
                #
                
                #
                # Next, add subgraphs by increasing distance from the seed
                # subgraph, until the maximum allowed number of nodes per region
                # graph is exceeded or until all nodes are in the decomposition
                #
                
                for i_subgraph_candidate_to_add in subgraph_indeces_sorted_by_distance_from_the_seed_subgraph[1:]:

                        #
                        # Only consider subgraphs that have not yet
                        # been added to the decomposition
                        #
                                
                        if i_subgraph_candidate_to_add in indeces_of_subgraphs_used_for_this_decomposition:

                                continue


                        the_subgraph_to_add = list_of_subgraphs[i_subgraph_candidate_to_add]
                        the_subgraph_to_add_nodes = list(the_subgraph_to_add.nodes)

                
                        #
                        # Remember to first subtract out the nodes
                        # that are already in the decomposition
                        #
                                
                        the_subgraph_to_add_nodes = list(set(the_subgraph_to_add_nodes).difference(set(nodes_in_decomposition)))                        

                
                        #
                        # Then, add the subgraph to the region graph
                        #
                        
                        this_region_graph.extend(the_subgraph_to_add_nodes)
                        indeces_of_subgraphs_used_for_this_decomposition.append(i_subgraph_candidate_to_add)

                        
                        #
                        # Keep track of which nodes are already in the
                        # decomposition
                        #

                        nodes_in_decomposition.extend(the_subgraph_to_add_nodes)

                        #
                        # Keep track of which nodes still need to be added
                        # to the decomposition
                        #
                                
                        nodes_remaining = list(set(nodes_remaining).difference(set(the_subgraph_to_add_nodes)))

                        #
                        # Stop adding subgraphs to this region graph when the region graph is ready, or when
                        # all nodes in the real-life graph are in the decomposition
                        #
                        
                        if len(this_region_graph) >= parameters.maximum_number_of_nodes_per_region_graph or len(nodes_remaining) == 0:

                                break

                #
                # Add the new region graph to the decomposition
                #

                decomposition.append(this_region_graph)


        #
        # The decomposition is now ready. Check that each node is
        # included in the final, decomposed graph exactly once.
        #

        all_nodes_in_decomposition = []
        
        for this_region_graph in decomposition:

                all_nodes_in_decomposition.extend(this_region_graph)

                
        unique_elements, counts = np.unique(all_nodes_in_decomposition, return_counts = True)

        if set(unique_elements) != set(list_of_all_nodes):

                print("ERROR! Not all nodes are included in the decomposition. Exiting.")
                exit(1)

        for c in counts:

                if c != 1:

                        print("ERROR! At least one node appeared a number of times not equal to one in the decomposition. Exiting.")
                        exit(1)
        

        print("Done.")

        
        #
        # Enable warnings again
        #
    
        warnings.resetwarnings()

        
        #
        # All done
        #
        
        return decomposition
        

#
# For a given path, compute the total weighted cost. Return this, as
# well as the pure length, duration, and ground damage components and
# their normalized values for the path. Return also the number of
# different assortments for each tour in the path.
#

def compute_cost_of_solution(this_solution, distance_matrix, ground_damage_matrix, wood_matrix, estimated_mean_length_of_path, estimated_mean_duration_of_path, estimated_mean_ground_damage_of_path):

        
        path = this_solution.path       
        length = 0.0
        duration = 0.0
        ground_damage = 0.0
        total_weighted_cost = 0.0

        
        #
        # Extract the number of different types of wood for each load
        # into this list
        #
        
        n_types_of_wood_by_load = []

        
        #
        # Iterate over the path, computing the cumulative cost as you
        # go
        #

        current_load = []
        
        for i in range(0, len(path) - 1):

                
                this_node = path[i]
                next_node = path[i+1]

                
                #
                # When departing from the unloading node, unload the
                # entire load and compute the duration of the
                # unloading
                #
                
                if this_node == 0 and len(current_load) > 0:

                        duration = duration + compute_duration_of_unloading(current_load)

                        #
                        # Also, gather data on the number of unique
                        # types of wood for each load
                        #
                        
                        this_n_types_of_wood = len(np.unique(np.array([pile[0] for pile in current_load])))
                        n_types_of_wood_by_load.append(this_n_types_of_wood)
                        
                        current_load = []

                
                #
                # When departing from a pickup node, load up the contents of that node
                #
                
                if this_node != 0:
                        
                        wood_type_to_load = wood_matrix[this_node, 0]
                        wood_amount_to_load = wood_matrix[this_node, 1]
                        current_load.append([wood_type_to_load, wood_amount_to_load])
                
                
                #
                # Add length, duration, and ground damage for
                # traversing the edge from this node to the next
                #

                length = length + distance_matrix[this_node, next_node]

                duration = duration + distance_matrix[this_node, next_node] / parameters.forwarder_speed
                
                current_load_volume = np.sum(np.array([pile[1] for pile in current_load]))
                ground_damage = ground_damage + ground_damage_matrix[this_node, next_node]*(1.0 + current_load_volume / parameters.ground_damage_model_v0)
                
                
                #
                # Before returning to the unloading node, compute
                # duration of loading the entire load so far
                #
                
                if next_node == 0:

                        duration = duration + compute_duration_of_loading(current_load)

                        
        
        #
        # Add duration for unloading the final load
        #

        duration = duration + compute_duration_of_unloading(current_load)

        
        #
        # Also, get the number of unique types of wood for the final
        # load
        #
        
        this_n_types_of_wood = len(np.unique(np.array([pile[0] for pile in current_load])))
        n_types_of_wood_by_load.append(this_n_types_of_wood)

        
        #
        # Compute total cost of the path using the given weights
        #

        #
        # Normalize the components first
        #

        standardized_length_component = length / estimated_mean_length_of_path
        standardized_duration_component = duration / estimated_mean_duration_of_path
        standardized_ground_damage_component = ground_damage / estimated_mean_ground_damage_of_path
        
        total_weighted_cost = parameters.w_length*standardized_length_component + parameters.w_duration*standardized_duration_component + parameters.w_ground_damage*standardized_ground_damage_component                

        #
        # All done
        #
        
        return total_weighted_cost, length, duration, ground_damage, standardized_length_component, standardized_duration_component, standardized_ground_damage_component, n_types_of_wood_by_load



#
# Compute duration (in s) of loading a given complete load
#

def compute_duration_of_loading(current_load):

        #
        # Get unique number of assortments to load
        #

        n_types_of_wood = len(np.unique(np.array([pile[0] for pile in current_load])))

        #
        # Use Manner (2013) Eq. (3) to get the loading time in seconds
        # per m**3
        #
        
        cost = parameters.t_load_a * np.power(n_types_of_wood, 2) + parameters.t_load_b * n_types_of_wood + parameters.t_load_c
        cost = cost*60.0

        #
        # Compute total volume of the load in m**3
        #
        
        total_volume_of_load = 0.0
        
        for pile in current_load:
                
                total_volume_of_load = total_volume_of_load + pile[1]

        #
        # Compute the total time required for the loading
        #
                
        cost = cost*total_volume_of_load
                
        return cost



#
# Compute duration (in s) of unloading a given complete load
#

def compute_duration_of_unloading(current_load):

        #
        # Get unique number of assortments to unload
        #

        n_types_of_wood = len(np.unique(np.array([pile[0] for pile in current_load])))

        #
        # Use Manner (2013) Eq. (4) to get the unloading time in
        # seconds per m**3
        #
        
        cost = parameters.t_unload_a * np.power(n_types_of_wood, 2) + parameters.t_unload_b * n_types_of_wood + parameters.t_unload_c
        cost = cost*60.0

        #
        # Compute total volume of the load in m**3
        #
        
        total_volume_of_load = 0.0
        
        for pile in current_load:
                
                total_volume_of_load = total_volume_of_load + pile[1]

        #
        # Compute the total time required for the unloading
        #
                
        cost = cost*total_volume_of_load
        
        return cost



#
# Construct the real-life graph of loading nodes, the unloading node,
# and the road segments between them.
#
# Each node is identified by its index, an integer running from 0 (the
# unloading or "roadside" node) up to N (the number of loading or
# "pickup" nodes).
#
# Each node is affiliated with an object:
#
# - smart_log_classes_aco.smart_log_pickup_node for loading nodes
#
# - smart_log_classes_aco.smart_log_roadside_node for the unloading
#   node
#
# The pickup node object holds the position (x, y) as well as the wood
# type and pile size (m**3) for that node. The roadside node object
# only holds the (x, y) position.
#
# The function does the following.
#
# First, read in a set of piles in the following format:
#
# <northing (m)> <easting (m)> <type of wood> <amount of wood (m**3)> <number of logs in the pile> <maximum log to log distance in the pile (m)>
#
# and read in a set of positions for the forwarder in the following
# format:
# 
# <altitude (m)> <northing (m)> <easting (m)>
#
# Then, create a road network from the given data, and project each
# pile to its nearest road segment.
#
# Next, from this information, construct a graph depicting the
# real-life locations of the piles and the assumed road segments
# between them.
#
# Also, plot the data at various stages of the process, if desired.
#
# Return the real-life graph and various lists giving the sizes and
# colors of nodes for plotting, as well as a dictionary of pile
# positions for plotting purposes.
#

def construct_real_life_graph_from_input_files(input_file_wood_piles, input_file_forwarder_position):

        #
        # Read in the data on wood piles
        #
        
        pile_data = np.loadtxt(input_file_wood_piles)

        #
        # Round the assortment types to the nearest integer just in case
        #
        
        pile_data[:, 2] = np.round(pile_data[:, 2])
        
        n_piles = pile_data.shape[0]

        list_of_wood_types = set([pile_data[i, 2].astype(int) for i in range(0, n_piles)])
                
        n_wood_types = len(list_of_wood_types)
        
        min_pile_size = np.min(pile_data[:, 3])
        max_pile_size = np.max(pile_data[:, 3])

        print("")
        print("Read in a total of %d lines of pile data of %d different assortments from file %s" % (n_piles, n_wood_types, input_file_wood_piles))
        print("")
        print("Number of logs in a pile ranges from %d to %d with a mean value of %f" % (np.min(pile_data[:, 4]), np.max(pile_data[:, 4]), np.mean(pile_data[:, 4])))
        print("")
        print("Pile volume ranges from %f to %f m**3 with a mean value of %f m**3" % (np.min(pile_data[:, 3]), np.max(pile_data[:, 3]), np.mean(pile_data[:, 3])))
        print("")
        print("Total volume of piles is %f m**3" % np.sum(pile_data[:, 3]))
        print("")
        print("The maximum separation of two logs in a pile runs from %f to %f m with a mean value of %f m" % (np.min(pile_data[:, 5]), np.max(pile_data[:, 5]), np.mean(pile_data[:, 5])))

        #
        # Find and print out the minimum and maximum easting and northing of the pile data
        #

        min_easting_pile_data = np.min(pile_data[:, 1])
        max_easting_pile_data = np.max(pile_data[:, 1])
        min_northing_pile_data = np.min(pile_data[:, 0])
        max_northing_pile_data = np.max(pile_data[:, 0])

        print("")
        print("The pile data runs from easting %f m to %f m and northing %f m to %f m" % (min_easting_pile_data, max_easting_pile_data, min_northing_pile_data, max_northing_pile_data))
        
        
        #
        # Create a dictionary for the wood types, i.e., the
        # assortments, where the original wood types in the data
        #
        # (142, 34, 21, ...)
        #
        # are remapped to a running index
        #
        # (1, 2, ...)
        #
        # The dictionary is thus of the format
        #
        # {original wood type : new wood type}.
        #
        
        wood_type_dictionary = {}
        i_wood_type = 1

        for wood_type in list_of_wood_types:
                
                wood_type_dictionary[wood_type] = i_wood_type
                i_wood_type = i_wood_type + 1


        #
        # Read in the data on forwarder position
        #

        forwarder_position_data = np.loadtxt(input_file_forwarder_position)
        n_positions = forwarder_position_data.shape[0]

        print("")
        print("Read in a total of %d lines of forwarder position data from file %s" % (n_positions, input_file_forwarder_position))

        #
        # These will be used for axes limits in plots as well as
        # limits for iterating windows over the data
        #
        
        min_x_forwarder_position_data = np.amin(forwarder_position_data[:, 2])
        max_x_forwarder_position_data = np.amax(forwarder_position_data[:, 2])
        min_y_forwarder_position_data = np.amin(forwarder_position_data[:, 1])
        max_y_forwarder_position_data = np.amax(forwarder_position_data[:, 1])
        
        #
        # Construct the unique road network, i.e., create a unique
        # line from the forwarder position data. Do this through the
        # following procedure:
        #
        # - Loop a window over the area
        #
        # - If you find enough points in a window, compute the average
        #   position of the points and set a road point there
        #

        print("")
        print("Now creating road from forwarder positional data...")
        
        unambiguous_road_points = []

        #
        # Loop over the window position in two dimensions, and find
        # the points denoting the road
        #
        
        for window_x in min_x_forwarder_position_data + np.arange(0, math.floor((max_x_forwarder_position_data - min_x_forwarder_position_data) / parameters.road_point_identification_window_size) + 1 + 1)*parameters.road_point_identification_window_size:
        
                for window_y in min_y_forwarder_position_data + np.arange(0, math.floor((max_y_forwarder_position_data - min_y_forwarder_position_data) / parameters.road_point_identification_window_size) + 1 + 1)*parameters.road_point_identification_window_size:
        
                        this_window = [window_x, window_x + parameters.road_point_identification_window_size, window_y, window_y + parameters.road_point_identification_window_size]

                        #
                        # Find points in this window
                        #
                        
                        these_points = forwarder_position_data[np.where(forwarder_position_data[:, 2] >= this_window[0])]
                        these_points = these_points[np.where(these_points[:, 2] < this_window[1])]
                        these_points = these_points[np.where(these_points[:, 1] >= this_window[2])]
                        these_points = these_points[np.where(these_points[:, 1] < this_window[3])]

                        #
                        # Cut out the altitude component
                        #
                        
                        these_points = these_points[:, 1:3]

                        #
                        # Switch to the regular xy-view
                        #
                        
                        these_points_temp = np.copy(these_points)
                        these_points[:, 0] = these_points_temp[:, 1]
                        these_points[:, 1] = these_points_temp[:, 0]

                        #
                        # If we have enough points in this window,
                        # make the mean position of these points a
                        # road point
                        #
                        
                        if these_points.shape[0] >= parameters.road_point_identification_min_n_points_per_window:

                                these_points_mean = np.mean(these_points, axis = 0)
                                unambiguous_road_points.append(these_points_mean)

        unambiguous_road_points = np.array(unambiguous_road_points)
        n_unambiguous_road_points = unambiguous_road_points.shape[0]

        #
        # Make sure that more than zero road points were created.
        #
        
        if n_unambiguous_road_points == 0:

                print("")
                print("ERROR! Created zero road points. Try and adjust the parameters for road creation.")
                print("Exiting.")
                print("")
                exit(1)

        print("Done. Created a total of %d road points." % n_unambiguous_road_points)

        
        #
        # Project the piles onto the road. Do this using the following
        # approach:
        #
        # - For each pile position, find the closest road point. We
        #   call this the primary road point.
        #
        # - Then, find the closest road points to this primary road
        #   point. We call these the secondary road points.
        #
        # - Fit a natural spline through the primary and secondary
        #   road points
        #
        # - Set the pile position to the point on the spline
        #   interpolation that is closest the original pile position
        #

        print("")
        print("Now projecting pile positions onto the road network...")

        #
        # Loop over the piles
        #

        new_pile_positions = np.zeros([pile_data.shape[0], 2])
        
        for i_pile in range(0, pile_data.shape[0]):

                #
                # Find the nearest road point to this pile, i.e., the
                # primary road point
                #
                
                this_pile_position = np.array([pile_data[i_pile, 1], pile_data[i_pile, 0]])
                
                road_points_and_distances = np.zeros([unambiguous_road_points.shape[0], 3])

                for j_pile in range(0, unambiguous_road_points.shape[0]):

                        this_road_point_x = unambiguous_road_points[j_pile, 0]
                        this_road_point_y = unambiguous_road_points[j_pile, 1]
                        this_distance = np.sqrt( (this_road_point_x - this_pile_position[0])**2 + (this_road_point_y - this_pile_position[1])**2 )
                        road_points_and_distances[j_pile, :] = np.array([this_road_point_x, this_road_point_y, this_distance])

                #
                # Sort the road points by distance and get the primary
                # road point
                #
                
                sorted_road_points_and_distances = road_points_and_distances[road_points_and_distances[:, 2].argsort()]
                primary_road_point = np.array([sorted_road_points_and_distances[0, 0], sorted_road_points_and_distances[0, 1]])

                #
                # Then, find the secondary road points
                #
                
                secondary_road_points_and_distances = np.zeros([unambiguous_road_points.shape[0], 3])

                for j_pile in range(0, unambiguous_road_points.shape[0]):
                        
                        this_road_point_x = unambiguous_road_points[j_pile, 0]
                        this_road_point_y = unambiguous_road_points[j_pile, 1]
                        this_distance = np.sqrt( (this_road_point_x - primary_road_point[0])**2 + (this_road_point_y - primary_road_point[1])**2 )
                        secondary_road_points_and_distances[j_pile, :] = np.array([this_road_point_x, this_road_point_y, this_distance])


                sorted_secondary_road_points_and_distances = secondary_road_points_and_distances[secondary_road_points_and_distances[:, 2].argsort()]
                secondary_road_points = sorted_secondary_road_points_and_distances[1:(parameters.n_road_points_for_pile_placement + 1), 0:2]

                #
                # Gather up the road points for which you will do the
                # spline interpolation, then interpolate.
                #

                road_points_for_interpolation = np.vstack((primary_road_point, secondary_road_points))
                
                #
                # Do the spline interpolation on the collected road
                # points. Start from the point closest to any corner
                # of the rectangle which just encompasses the road
                # points. Sort the points by distance from this corner
                # point, and define a parametric curve for which to do
                # the spline interpolation.
                #
                
                #
                # Find the road point closest to a corner of the
                # minimal rectangle just encompassing all the road
                # points
                #
                
                corner_points = np.zeros([4, 2])
                corner_points[0, :] = np.array([np.max(road_points_for_interpolation[:, 0]), np.max(road_points_for_interpolation[:, 1])]) # North-east corner
                corner_points[1, :] = np.array([np.min(road_points_for_interpolation[:, 0]), np.max(road_points_for_interpolation[:, 1])]) # North-west corner
                corner_points[2, :] = np.array([np.min(road_points_for_interpolation[:, 0]), np.min(road_points_for_interpolation[:, 1])]) # South-west corner
                corner_points[3, :] = np.array([np.max(road_points_for_interpolation[:, 0]), np.min(road_points_for_interpolation[:, 1])]) # South-east corner

                smallest_distance = np.Inf
                i_smallest_distance = 0
                
                for i_road_point in range(0, road_points_for_interpolation.shape[0]):
                        
                        for i_corner_point in range(0, corner_points.shape[0]):
                                
                                this_distance = np.linalg.norm(corner_points[i_corner_point, :] - road_points_for_interpolation[i_road_point, :])

                                if this_distance < smallest_distance:

                                        smallest_distance = this_distance
                                        i_smallest_distance = i_road_point
                
                road_point_closest_to_corner = np.array([road_points_for_interpolation[i_smallest_distance, 0], road_points_for_interpolation[i_smallest_distance, 1]])
                
                #
                # Sort the road points by distance from the corner point
                #

                road_points_for_interpolation_with_distances_from_corner_point = np.zeros([road_points_for_interpolation.shape[0], 3])
                
                for i_road_point in range(0, road_points_for_interpolation.shape[0]):

                        this_distance = np.linalg.norm(road_point_closest_to_corner - road_points_for_interpolation[i_road_point, :])
                        road_points_for_interpolation_with_distances_from_corner_point[i_road_point, :] = np.array([road_points_for_interpolation[i_road_point, 0], road_points_for_interpolation[i_road_point, 1], this_distance])
                
                road_curve_points = road_points_for_interpolation_with_distances_from_corner_point[road_points_for_interpolation_with_distances_from_corner_point[:, 2].argsort()]

                #
                # Use a generic parameter t value to parametrize the road curve
                #
                
                road_curve_points[:, 2] = np.arange(0, road_curve_points.shape[0])

                #
                # Fit a natural cubic spline to the parametric curve data (a, x(a), y(a))
                #
                
                road_spline = CubicSpline(road_curve_points[:, 2], road_curve_points[:, 0:2], bc_type = 'natural')

                #
                # Set the new pile position to the closest point on the spline
                #
                
                min_t = road_curve_points[0, 2]
                max_t = road_curve_points[-1, 2]
                t_points = np.linspace(min_t, max_t, (max_t - min_t) / parameters.road_spline_parameter_step)                
                spline_points = road_spline(t_points)
                
                i_closest_road_point_to_pile = np.argmin(np.linalg.norm(this_pile_position - spline_points, axis = 1))
                this_new_pile_position = spline_points[i_closest_road_point_to_pile, :]

                #
                # Save the new pile position
                #
                
                new_pile_positions[i_pile, :] = this_new_pile_position

                #
                # Plot
                #
                #plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                #plt.plot(road_points_for_interpolation[:, 0], road_points_for_interpolation[:, 1], 'rx', linewidth=0, markersize=10)
                #plt.plot(spline_points[:, 0], spline_points[:, 1], 'b-', linewidth=3)
                #plt.plot(this_pile_position[0], this_pile_position[1], 'mo', markersize=10)
                #plt.plot(this_new_pile_position[0], this_new_pile_position[1], 'ms', markersize=10)
                #plt.axis('equal')
                #plt.show()
                

        print("Done.")

        #
        # Compute statistics on the shifts made to pile positions upon
        # projecting them to the road network
        #

        old_pile_positions = np.vstack((pile_data[:, 1], pile_data[:, 0])).T
        pile_shifts = np.linalg.norm(new_pile_positions - old_pile_positions, axis = 1)
        
        min_shift_in_pile_positions = np.min(pile_shifts)
        mean_shift_in_pile_positions = np.mean(pile_shifts)
        max_shift_in_pile_positions = np.max(pile_shifts)
        std_of_shift_in_pile_positions = np.std(pile_shifts)
        
        print("")
        print("Shifts in pile positions upon projecting onto the road network ranged from %f m to %f m, with a mean of %f m and a standard deviation of %f m" % (min_shift_in_pile_positions, max_shift_in_pile_positions, mean_shift_in_pile_positions, std_of_shift_in_pile_positions))
        
        
        #
        # Construct the real-life graph consisting of the piles and
        # the connections between them
        #
        
        print("")
        print("Now constructing the full graph of the pile nodes...")

        #
        # 1. Add all the nodes to the graph. No connections yet.
        #

        #
        # This list holds the pick-up node objects that constitute the
        # data for each individual node
        #
        
        list_of_pickup_nodes = []
        
        for i_pile in range(0, n_piles):

                this_x_position = new_pile_positions[i_pile, 0]
                this_y_position = new_pile_positions[i_pile, 1]

                #
                # Wood type in the graph object is set here to the new coding, i.e., 1, 2, ...
                #
                
                this_wood_type = wood_type_dictionary[int(pile_data[i_pile, 2])]

                this_pile_size = pile_data[i_pile, 3]
	
                list_of_pickup_nodes.append(smart_log_classes_aco.smart_log_pickup_node(this_x_position, this_y_position, this_wood_type, this_pile_size))

        #
        # Create an empty graph
        #
        
        the_graph = networkx.Graph()

        #
        # Add the roadside node. Label it 0.
        #
        
        the_graph.add_node(0, data = smart_log_classes_aco.smart_log_roadside_node(parameters.x_roadside_position, parameters.y_roadside_position))

        #
        # Add the pickup nodes. Label them as 1, 2, ...
        #
        
        i_node = 1
        
        for node_data in list_of_pickup_nodes:

                the_graph.add_node(i_node, data = node_data)

                i_node = i_node + 1

                
        #
        # 2. Then, add the connections, i.e., the edges. The idea is
        # to form a string of pick-up nodes, which approximates the
        # road network and the piles along that road network. Do this
        # by the following approach:
        #
        # - For each node, find its nearest neighbor and connect to it
        #
        # - Then, find the nearest neighbor in a narrow sector
        #   directly opposite to the first nearest neighbor
        #
        #
        # TODO: Parallellize this?
        #

        #
        # First, add an edge from the roadside node simply to the
        # pick-up node closest to it
        #
        
        roadside_node_position = np.array([the_graph.nodes[0]['data'].x_position, the_graph.nodes[0]['data'].y_position])
        
        smallest_distance_so_far = np.Inf
        i_nearest_pile = -1
        
        for i_pile in range(0, n_piles):

                this_pile_position = np.array([new_pile_positions[i_pile, 0], new_pile_positions[i_pile, 1]])
                this_distance = np.linalg.norm(this_pile_position - roadside_node_position)

                if this_distance < smallest_distance_so_far:
                        
                        smallest_distance_so_far = this_distance
                        i_nearest_pile = i_pile

        the_graph.add_edge(0, i_nearest_pile + 1)

        #
        # BEGIN DEBUG
        #

        #print("")
        #print("The nearest pile to the roadside point was %d " % (i_nearest_pile + 1))

        #
        # END DEBUG
        #
        
        #
        # Then, add the edges between the pick-up nodes
        #
        
        this_pile_neighbor_positions = None
        
        for i_pile in range(0, n_piles):

                #
                # Skip the pile nearest to the roadside node, as it
                # already has the one connection it needs
                #
                
                if i_pile == i_nearest_pile:

                        continue

                this_pile_position = np.array([new_pile_positions[i_pile, 0], new_pile_positions[i_pile, 1]])

                #
                # BEGIN DEBUG
                #

                #print("")
                #print("Now handling node number %d" % (i_pile + 1))
                
                #
                # END DEBUG
                #
                
                
                #
                # Find the nearest neighbor for this pile
                #
                
                smallest_distance_so_far = np.Inf
                i_nn = -1
                
                for j_pile in range(0, n_piles):

                        if j_pile != i_pile:

                                this_neighboring_pile_position = np.array([new_pile_positions[j_pile, 0], new_pile_positions[j_pile, 1]])
                                this_distance = np.linalg.norm(this_pile_position - this_neighboring_pile_position)

                                #
                                # The potential neighbor must be within an allowed distance range
                                #
                                
                                if this_distance > parameters.pile_connectivity_max_distance_between_neighbors:
                                        
                                        continue
                                
                                if this_distance < smallest_distance_so_far:
                                        
                                        smallest_distance_so_far = this_distance
                                        i_nn = j_pile


                if i_nn < 0:
                        
                        print("ERROR! No nearest neighbor found for pile %d. Exiting." % (i_pile + 1))
                        exit(1)
                        
                                        
                this_pile_nearest_neighbor_position = np.array([new_pile_positions[i_nn, 0], new_pile_positions[i_nn, 1]])
                pile_to_nearest_neighbor = this_pile_nearest_neighbor_position - this_pile_position
                
                #
                # Add edge from this node to its nearest neighbor
                #
                
                the_graph.add_edge(i_pile + 1, i_nn + 1)

                #
                # BEGIN DEBUG
                #

                #print("Found nearest neighbor %d" % (i_nn + 1))

                #
                # END DEBUG
                #
                
                
                #
                # Then, find the nearest neighbor in the allowed
                # sector. We call this the second-nearest neighbor. If
                # the second-nearest neighbor is not found, increase
                # the sector. Repeat until you reach the maximum
                # allowed sector size.
                #
                
                i_second_nn = -1
                sector_range = parameters.pile_connectivity_initial_sector_range

                while i_second_nn < 0:
                        
                        smallest_distance_so_far = np.Inf
                
                        for j_pile in range(0, n_piles):

                                if j_pile != i_pile and j_pile != i_nn:

                                        this_neighboring_pile_position = np.array([new_pile_positions[j_pile, 0], new_pile_positions[j_pile, 1]])
                                        this_distance = np.linalg.norm(this_pile_position - this_neighboring_pile_position)

                                        #
                                        # The potential neighbor must be within an allowed distance range
                                        #
                                        
                                        if this_distance > parameters.pile_connectivity_max_distance_between_neighbors:
                                                
                                                continue
                                        
                                        #
                                        # Compute the angle to the potential neighbor
                                        #
                                        
                                        pile_to_second_nearest_neighbor = this_neighboring_pile_position - this_pile_position

                                        arccos_input_nominator = np.dot(pile_to_nearest_neighbor, pile_to_second_nearest_neighbor)
                                        arccos_input_denominator = np.linalg.norm(pile_to_nearest_neighbor) * np.linalg.norm(pile_to_second_nearest_neighbor)

                                        #
                                        # If the denominator is zero, just find the nearest neighbor at any angle
                                        #
                                        
                                        if arccos_input_denominator == 0:

                                                if this_distance < smallest_distance_so_far:
                                                
                                                        smallest_distance_so_far = this_distance
                                                        i_second_nn = j_pile

                                                continue

                                        #
                                        # Otherwise, see if this potential neighbor is in the allowed sector
                                        #
                                        
                                        arccos_input = arccos_input_nominator / arccos_input_denominator

                                        #
                                        # Clip the input for the arccos function, as numerical noise may cause a value higher than 1.0 in absolute value
                                        #
                                        
                                        arccos_input = np.clip(arccos_input, -1.0, 1.0)
                                        this_angle = np.arccos(arccos_input)

                                        if (np.pi - this_angle) * 180.0 / np.pi < sector_range:
                                
                                                if this_distance < smallest_distance_so_far:
                                                
                                                        smallest_distance_so_far = this_distance
                                                        i_second_nn = j_pile
                                                
                        if i_second_nn < 0:

                                sector_range = sector_range + parameters.pile_connectivity_sector_range_increase

                                if sector_range > parameters.pile_connectivity_max_sector:

                                        print("Warning: No second-nearest neighbor found for pile %d at (%f, %f)" % (i_pile + 1, this_pile_position[0], this_pile_position[1]))
                                        this_pile_neighbor_positions = np.array([this_pile_position, this_pile_nearest_neighbor_position])
                                        break
                                
                        else:

                                this_pile_second_nearest_neighbor_position = np.array([new_pile_positions[i_second_nn, 0], new_pile_positions[i_second_nn, 1]])
                                this_pile_neighbor_positions = np.array([this_pile_second_nearest_neighbor_position, this_pile_position, this_pile_nearest_neighbor_position])

                                #
                                # Add edge from this node to its second-nearest neighbor
                                #

                                the_graph.add_edge(i_pile + 1, i_second_nn + 1)

                                #
                                # BEGIN DEBUG
                                #
                                
                                #print("Found second-nearest neighbor %d" % (i_second_nn + 1))

                                #
                                # END DEBUG
                                #
                                
        #
        # 3. Now that you've added the nodes and edges to the graph,
        # compute the edge lengths, i.e., the eucledian length of each
        # edge.
        #

        for e in the_graph.edges.items():

                this_edge_begin_node = e[0][0]
                this_edge_end_node = e[0][1]
                
                x_begin = the_graph.nodes[this_edge_begin_node]['data'].x_position
                y_begin = the_graph.nodes[this_edge_begin_node]['data'].y_position
                x_end = the_graph.nodes[this_edge_end_node]['data'].x_position
                y_end = the_graph.nodes[this_edge_end_node]['data'].y_position

                this_edge_length = np.sqrt(np.power(x_end-x_begin, 2) + np.power(y_end-y_begin, 2))

                e[1]['length'] = this_edge_length
        
        #
        # The real-life graph is now ready
        #
                
        print("Done. Created a real-life graph with a total of %d nodes and %d edges." % (len(list(the_graph.nodes)), len(list(the_graph.edges))))

        
        #
        # Visualize the results
        #

        #
        # Create a dictionary of colors of the format
        #
        # {original wood type : color}
        #
        # This mapping of wood type to color will be used throughout
        # the program.
        #
        
        wood_type_color_space = cm.rainbow(np.linspace(0, 1, n_wood_types))
        
        wood_type_color_dictionary = {}
        
        for wood_type in wood_type_dictionary:
                
                wood_type_color_dictionary[wood_type] = wood_type_color_space[wood_type_dictionary[int(wood_type)] - 1]
        
        
        #
        # Create a list of pile colors that matches the pile data
        # element by element
        #
        
        pile_marker_colors = []
        
        for wood_type in pile_data[:, 2].astype(int):

                this_pile_marker_color = wood_type_color_dictionary[wood_type]
                pile_marker_colors.append(this_pile_marker_color)

        
        #
        # Create a list of pile sizes that matches the pile data
        # element by element
        #

        pile_marker_sizes = []
        
        for pile_size in pile_data[:, 3]:

                if(max_pile_size - min_pile_size > 0):

                        this_pile_marker_size = (((pile_size - min_pile_size) / (max_pile_size - min_pile_size) + 1.0) * parameters.pickup_node_marker_scaling_factor)**3

                else:

                        this_pile_marker_size = parameters.default_pile_marker_size
                        
                pile_marker_sizes.append(this_pile_marker_size)

        
        # DEBUG
        #xmin = min_x_forwarder_position_data + np.random.random()*(max_x_forwarder_position_data - min_x_forwarder_position_data)
        #xmax = xmin + 50.0
        #ymin = min_y_forwarder_position_data + np.random.random()*(max_y_forwarder_position_data - min_y_forwarder_position_data)
        #ymax = ymin + 50.0
        #these_axis_lims = [xmin, xmax, ymin, ymax]
        #
        
        # Set a reasonable default font size for all plots
        plt.rc('font', size = 15)


        #
        # Create a pie chart of wood type by number of logs in each
        # assortment
        #

        wood_types, wood_type_counts = np.unique(pile_data[:, 2].astype(int), return_counts = True)
        
        pie_labels = []
        pie_wedge_colors = []

        for wt in wood_types:

                pie_labels.append('Assortment ' + str(int(wt)))
                pie_wedge_colors.append(wood_type_color_dictionary[int(wt)])
                
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.pie(wood_type_counts, labels = pie_labels, colors = pie_wedge_colors, autopct = '%1.1f%%')
        plt.title('Distribution of assortments by number of piles', pad = 25)

        
        if parameters.save_figures_to_files:
        
                plt.savefig('pie_chart_of_wood_type_distribution.png')

        plt.close()



        #
        # Create a histogram of pile volume
        #
    
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = 15)
        plt.rc('axes', titlesize = 25)

        pile_sizes = pile_data[:, 3]
    
        #
        # Set the bins for the histogram. The left edge of the first bin is at zero.
        #

        bins = 0.0 + np.arange(0, math.floor((np.max(pile_sizes) + 2.0*parameters.pile_volume_histogram_bin_step - 0.0) / parameters.pile_volume_histogram_bin_step) + 1)*parameters.pile_volume_histogram_bin_step
    
        #
        # Then plot the histogram
        #
    
        plt.hist(pile_sizes, bins)
        plt.title('Distribution of pile size, total number of piles = ' + str(pile_sizes.shape[0]), pad = parameters.title_pad)
        plt.xlabel('Pile size (m$^3$)')
        plt.ylabel('Occurrencies')
    
        if parameters.save_figures_to_files:
        
                plt.savefig('histogram_of_pile_volume.png')

        plt.close()
        
    
        #
        # Create a histogram of number of logs in a pile
        #
    
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = 15)
        plt.rc('axes', titlesize = 25)

        pile_sizes = pile_data[:, 4]
    
        #
        # Set the bins for the histogram. The left edge of the first bin is at 0.5.
        #

        bins = 0.5 + np.arange(0, math.floor((np.max(pile_sizes) + 2.0*parameters.pile_number_histogram_bin_step - 0.0) / parameters.pile_number_histogram_bin_step) + 1)*parameters.pile_number_histogram_bin_step
    
        #
        # Then plot the histogram
        #
    
        plt.hist(pile_sizes, bins)
        plt.title('Distribution of pile size, total number of piles = ' + str(pile_sizes.shape[0]), pad = parameters.title_pad)
        plt.xlabel('Number of logs in pile')
        plt.ylabel('Occurrencies')
    
        if parameters.save_figures_to_files:
        
                plt.savefig('histogram_of_number_of_logs_in_a_pile.png')

        plt.close()

    
        #
        # Create a histogram of largest log to log distance in a pile
        #
    
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = 15)
        plt.rc('axes', titlesize = 25)

        pile_sizes = pile_data[:, 5]
    
        #
        # Set the bins for the histogram. The left edge of the first bin is at zero.
        #

        bins = 0.0 + np.arange(0, math.floor((np.max(pile_sizes) + 2.0*parameters.pile_extent_histogram_bin_step - 0.0) / parameters.pile_extent_histogram_bin_step) + 1)*parameters.pile_extent_histogram_bin_step
    
        #
        # Then plot the histogram
        #
    
        plt.hist(pile_sizes, bins)
        plt.title('Distribution of pile size, total number of piles = ' + str(pile_sizes.shape[0]), pad = parameters.title_pad)
        plt.xlabel('Maximum distance of logs in pile (m)')
        plt.ylabel('Occurrencies')
	
        if parameters.save_figures_to_files:
            
                plt.savefig('histogram_of_pile_extent.png')

        plt.close()
        

        #
        # Create a heat map of the pile locations for each wood type,
        # i.e., assortment
        #

        bin_edges_x = np.array(np.linspace(min_x_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_x_forwarder_position_data + parameters.margin_for_pile_and_road_plots, parameters.n_bins_wood_type_histogram))
        bin_edges_y = np.array(np.linspace(min_y_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_y_forwarder_position_data + parameters.margin_for_pile_and_road_plots, parameters.n_bins_wood_type_histogram))

        #
        # First, get the minimum and maximum counts of wood locations
        # over all bins
        #

        min_counts = np.Inf
        max_counts = -1

        for wt in wood_type_dictionary:
                
                these_wood_type_locations = new_pile_positions[np.where(pile_data[:, 2].astype(int) == wt)]                        

                this_histogram_data, ex, ey = np.histogram2d(these_wood_type_locations[:, 0], these_wood_type_locations[:, 1], bins = [bin_edges_x, bin_edges_y])
                this_min_counts = np.min(this_histogram_data)
                this_max_counts = np.max(this_histogram_data)

                if this_min_counts < min_counts:

                        min_counts = this_min_counts

                if this_max_counts > max_counts:

                        max_counts = this_max_counts

        #
        # Then, plot the 2D histogram of wood locations separately for
        # each wood type
        #

        n_subplots_per_dimension = int(math.ceil(np.sqrt(n_wood_types)))

        #
        # Set a small font size for all text in the subplots to come
        #
        
        plt.rc('font', size = parameters.small_font_size)                

        fig, axes = plt.subplots(nrows = n_subplots_per_dimension, ncols = n_subplots_per_dimension, figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.title('Distribution of new locations of each wood type')

        n_subplots = 0

        #
        # Do this so the flattening operation succeeds even when
        # there's just a single wood type
        #
        
        if n_wood_types == 1:

                axes = np.array([axes])
        
        for wt, ax in zip(wood_type_dictionary, axes.flat):

                these_wood_type_locations = new_pile_positions[np.where(pile_data[:, 2].astype(int) == wt)]
                h, ex, ey, im  = ax.hist2d(these_wood_type_locations[:, 0], these_wood_type_locations[:, 1], bins = [bin_edges_x, bin_edges_y], vmin = min_counts, vmax = max_counts)
                ax.set_xlabel('Easting (m)')
                ax.set_ylabel('Northing (m)')
                ax.set_title('Assortment ' + str(wt) + '\n')
                ax.set_aspect('equal', adjustable = 'box')
        
                n_subplots = n_subplots + 1

        #
        # Revert to the default font size
        #
        
        plt.rc('font', size = parameters.default_font_size)

        #
        # Delete any remaining empty subplots
        #
        
        for i_ax in range(n_subplots_per_dimension**2, n_subplots, -1):
                axes.flat[i_ax-1].set_visible(False)

        fig.subplots_adjust(right = 0.8, hspace = 1.0, wspace = 0.4)
        cbar_axes = fig.add_axes([0.85, 0.15, 0.035, 0.6])
        cbar = fig.colorbar(im, cax = cbar_axes)
        cbar.set_label('Counts')


        if parameters.save_figures_to_files:
                
                plt.savefig('heat_map_of_piles_by_type.png')

        plt.close()
        
                
        #
        # Create a heat map of the pile locations regardless of wood
        # type, i.e., assortment
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.hist2d(new_pile_positions[:, 0], new_pile_positions[:, 1], bins = [bin_edges_x, bin_edges_y])
        cbar = plt.colorbar()
        cbar.set_label('Counts')
        plt.title('Distribution of new locations of wood')
        plt.xlabel('Easting (m)')        
        plt.ylabel('Northing (m)')
        plt.xlim([min_x_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_x_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.ylim([min_y_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_y_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:
                
                plt.savefig('heat_map_of_piles_regardless_of_type.png')

        plt.close()
        

        #
        # Plot the raw forwarder positional data and the road data
        # derived from this
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.title('Road data from forwarder position data')
        plt.xlabel('Easting (m)')        
        plt.ylabel('Northing (m)')
        plt.xlim([min_x_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_x_forwarder_position_data + parameters.margin_for_pile_and_road_plots])               
        plt.ylim([min_y_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_y_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.plot(forwarder_position_data[:, 2], forwarder_position_data[:, 1], 'bo', markersize = parameters.raw_forwarder_data_markersize, fillstyle = 'none', linewidth = 0, markeredgewidth = parameters.raw_forwarder_position_point_edge_width)
        plt.plot(unambiguous_road_points[:, 0], unambiguous_road_points[:, 1], 'rx', markersize = parameters.road_point_markersize, linewidth = 0, markeredgewidth = parameters.road_point_marker_edge_width)

        road_data_legend_handles = [mlines.Line2D([], [], color='r', marker = 'x', linestyle='None', markersize = parameters.road_point_markersize, markeredgewidth = parameters.road_point_marker_edge_width, label = 'Derived road points'), mlines.Line2D([], [], color='b', marker = 'o', fillstyle = 'none', linestyle='None', markersize = parameters.raw_forwarder_data_markersize, label = 'Raw forwarder positional data')]

        plt.legend(handles = road_data_legend_handles)
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:
                
                plt.savefig('road_data.png')


        plt.close()
                

        #
        # Plot the old pile positions with the raw forwarder position
        # data
        #
        
        #
        # First, create a solid list of legend handles that matches
        # each pile marker to its correct color. This will be useful
        # for visualizations throughout the code.
        #
        
        pile_color_legend_handles = []
        
        for wood_type in wood_type_dictionary:
                
                this_pile_marker_color = wood_type_color_dictionary[wood_type]
                
                pile_color_legend_handles.append(mlines.Line2D([], [], color = this_pile_marker_color, marker = '.', linestyle = 'None', markersize = parameters.legend_wood_pile_marker_size, label = 'Assortment ' + str(wood_type)))
        
        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))

        # DEBUG
        #plt.axis('equal')
        plt.title('Pile positions and road network (raw data)')
        plt.xlabel('Easting (m)')        
        plt.ylabel('Northing (m)')
        plt.xlim([min_x_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_x_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.ylim([min_y_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_y_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        #plt.xlim(these_axis_lims[0:2])
        #plt.ylim(these_axis_lims[2:4])
        plt.scatter(pile_data[:, 1], pile_data[:, 0], color = pile_marker_colors, sizes = pile_marker_sizes)
        plt.plot(forwarder_position_data[:, 2], forwarder_position_data[:, 1], 'k-', linewidth = 1)
        # DEBUG
        #plt.plot(unambiguous_road_points[:, 0], unambiguous_road_points[:, 1], 'rx', markersize = parameters.road_point_markersize, linewidth = 0, markeredgewidth = parameters.road_point_marker_edge_width)
        plt.legend(handles = pile_color_legend_handles)
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:

                plt.savefig('raw_forwarder_positional_data_with_original_pile_positions.png')

        plt.close()
        

        #
        # Plot the new pile positions with the raw forwarder position
        # data
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        # DEBUG
        #plt.axis('equal')
        plt.title('Pile positions and road network (piles projected onto road)')
        plt.xlabel('Easting (m)')
        plt.ylabel('Northing (m)')
        #plt.xlim(these_axis_lims[0:2])
        #plt.ylim(these_axis_lims[2:4])
        plt.xlim([min_x_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_x_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.ylim([min_y_forwarder_position_data - parameters.margin_for_pile_and_road_plots, max_y_forwarder_position_data + parameters.margin_for_pile_and_road_plots])
        plt.scatter(new_pile_positions[:, 0], new_pile_positions[:, 1], color = pile_marker_colors, sizes = pile_marker_sizes)
        plt.plot(forwarder_position_data[:, 2], forwarder_position_data[:, 1], 'k-', linewidth = 1, markersize = 5)
        # DEBUG
        #plt.plot(unambiguous_road_points[:, 0], unambiguous_road_points[:, 1], 'rx', markersize = parameters.road_point_markersize, linewidth = 0, markeredgewidth = parameters.road_point_marker_edge_width)
        plt.legend(handles = pile_color_legend_handles)
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:

                plt.savefig('raw_forwarder_positional_data_with_new_pile_positions.png')

        plt.close()
        
                
        #
        # Finally, plot the created real-life graph
        #

        #
        # Create a dictionary for the positions of the nodes. This
        # will be useful for visualizations throughout the code.
        #
        
        positions = {}

        for n in the_graph.nodes.items():
                
                positions[n[0]] = [n[1]['data'].x_position, n[1]['data'].y_position]

        #
        # Create a list of node colors to use for plots of graphs
        # including the roadside node. Prepend a color for the
        # roadside node (light grey). This will be useful for
        # visualizations throughout the code.
        #

        pile_node_colors = pile_marker_colors.copy()
        pile_node_colors.insert(0, [0.82, 0.82, 0.82])

        #
        # Create a list of node sizes to use for plots of graphs
        # including the roadside node. Prepend a size for the roadside
        # node. This will be useful for visualizations throughout the
        # code.
        #

        pile_node_sizes = pile_marker_sizes.copy()
        pile_node_sizes.insert(0, parameters.roadside_node_scale)

        #
        # Plot
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = 15)
        plt.title('Real-life graph')

        
        #
        # Disable all warnings for this plot, to avoid unnecessary deprecation warnings
        #
        
        warnings.simplefilter('ignore')

        networkx.draw(the_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, with_labels = parameters.draw_graph_nodes_with_labels)

        #
        # Enable warnings again
        #
        
        warnings.resetwarnings()

        plt.legend(handles = pile_color_legend_handles)

        #
        # Get a set of axis limits to use for all graph plots from now on
        #

        xmin, xmax, ymin, ymax = plt.axis()
        graph_axis_limits = [xmin, xmax, ymin, ymax]
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:

                plt.savefig('real-life_graph.png')


        plt.close()

        
        #
        # Check that any node is accessible from any other node. If not, give error and exit.
        #

        print("")
        print("Now checking that each node in the graph is accessible from any other node...")

        if not networkx.is_connected(the_graph):

                print("ERROR! The graph is not connected. Exiting.")
                exit(1)
        
        print("Done.")


        #
        # Compute total length of the road network
        #

        total_length_of_road_network = 0.0
        
        for e in the_graph.edges.items():

                total_length_of_road_network = total_length_of_road_network + e[1]['length']

        
        print("")
        print("Total length of the road network is %f m" % total_length_of_road_network)


        #
        # Finally, return the created real-life graph along with
        # objects needed for creating consistent plots in other parts
        # of the code later on
        #
        
        return the_graph, pile_node_sizes, pile_node_colors, pile_color_legend_handles, wood_type_color_dictionary, wood_type_dictionary, graph_axis_limits, positions



#
# Compute and return the static cost matrices needed for the
# optimization. Do the computation for a given list of nodes
# in a connected graph. The cost matrices are:
#
# - The distance matrix D_ij
#
# - The ground damage matrix G_ij
#
# Return also the mean harvestability index over the computed set of
# nodes.
#

def compute_static_cost_matrices(the_graph, list_of_nodes, harvestability_map, harvestability_map_final_boundingbox_min_easting, harvestability_map_final_boundingbox_max_easting, harvestability_map_final_boundingbox_min_northing, harvestability_map_final_boundingbox_max_northing):

        #
        # Loop over all possible pairs of nodes (i, j) in the list of
        # nodes given as input. For each pair:
        #
        # - Find the shortest path between the two in the real-life
        #   graph, and set the distance matrix element D_ij to the
        #   length L of this shortest path.
        #
        # - Find the mean harvestability H over the nodes in this
        #   path, and set the ground damage matrix G_ij to H*L.
        #

        #
        # Initialize the cost matrices. Set the default value for the
        # cost to infinity.
        #
        
        total_number_of_nodes = len(list(the_graph.nodes))

        the_distance_matrix = np.ones((total_number_of_nodes, total_number_of_nodes))*np.Inf

        the_ground_damage_matrix = np.ones((total_number_of_nodes, total_number_of_nodes))*np.Inf

        mean_harvestability = 0.0

        
        #
        # Compute the distance and ground damage matrix elements by
        # iterating all pairs of nodes
        #
        
        all_possible_pairs_of_nodes = itertools.combinations(list_of_nodes, 2)
        
        for this_node_pair in all_possible_pairs_of_nodes:
                
                this_begin_node = this_node_pair[0]
                this_end_node = this_node_pair[1]

                #
                # Find the shortest path and path length between these
                # nodes
                #
                
                this_shortest_path = networkx.shortest_path(the_graph, this_begin_node, this_end_node, weight = 'length', method = 'dijkstra')
                this_shortest_path_length = networkx.shortest_path_length(the_graph, this_begin_node, this_end_node, weight = 'length', method = 'dijkstra')

                #
                # Some nodes may have ended up directly on top of each
                # other when projected onto the road. Set a minimum
                # distance slightly greater than zero for such nodes.
                #

                if this_shortest_path_length == 0.0 and this_begin_node != this_end_node:
                
                        print("Warning: Setting distance matrix element (%d, %d) to %f" % (this_begin_node, this_end_node, parameters.minimum_cost_in_distance_matrix))
                        this_shortest_path_length = parameters.minimum_cost_in_distance_matrix


                #
                # Compute the mean harvestability over the nodes in
                # the obtained shortest path
                #

                list_of_harvestability = []

                for n in this_shortest_path:

                        this_node_attributes = the_graph.nodes[n]

                        this_node_x_position = this_node_attributes['data'].x_position
                        this_node_y_position = this_node_attributes['data'].y_position

                        i_harvest = (np.round(-1.0*(this_node_y_position - harvestability_map_final_boundingbox_max_northing) / parameters.harvestability_map_pixel_size)).astype(int)
                        j_harvest = (np.round((this_node_x_position - harvestability_map_final_boundingbox_min_easting) / parameters.harvestability_map_pixel_size)).astype(int)
                        this_harvestability = harvestability_map[i_harvest, j_harvest]

                        list_of_harvestability.append(this_harvestability)

                this_shortest_path_harvestability = np.mean(np.array(list_of_harvestability))
                
                #
                # Store the computed values for the matrix elements
                #
                        
                the_distance_matrix[this_begin_node, this_end_node] = this_shortest_path_length
                the_distance_matrix[this_end_node, this_begin_node] = this_shortest_path_length

                the_ground_damage_matrix[this_begin_node, this_end_node] = this_shortest_path_harvestability*this_shortest_path_length
                the_ground_damage_matrix[this_end_node, this_begin_node] = this_shortest_path_harvestability*this_shortest_path_length

                # DEBUG
                # print("Shortest path between nodes %d and %d which has a length of %f and harvestability of %f: " % (this_begin_node, this_end_node, this_shortest_path_length, this_shortest_path_harvestability), this_shortest_path)
                # END DEBUG
                
                continue
        
                ### BEGIN DEBUG

                print("Shortest path between nodes %d and %d which has a length of %f and harvestability of %f: " % (this_begin_node, this_end_node, this_shortest_path_length, this_shortest_path_harvestability), this_shortest_path)

                #
	        # Create a dictionary for the positions of the nodes
	        #
                positions = {}
        
                for n in the_graph.nodes.items():
                        
                        positions[n[0]] = [n[1]['data'].x_position, n[1]['data'].y_position]


                # Use a directed graph so that we can draw arrows
                diGraphForPath = networkx.DiGraph()
                diGraphForPath.add_nodes_from(the_graph)

                # Get the path edges for plotting
                this_shortest_path_edges = [(this_shortest_path[i], this_shortest_path[i+1]) for i in range (0, len(this_shortest_path) - 1)]
                
                #
                # Plot
                #

                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                plt.rc('font', size = parameters.default_font_size)
                plt.title('Real-life graph with shortest path between nodes %d and %d' % (this_begin_node, this_end_node))

                #
                # Draw just the nodes first. Draw the beginning and end nodes of the path at a larger scale and different color than the rest of the nodes.
                #                

                node_colors = ['g' for i_node in range(0, len(list(the_graph.nodes)))]
                node_colors[this_begin_node] = 'b'
                node_colors[this_end_node] = 'r'

                node_sizes = [0.0 for i_node in range(0, len(list(the_graph.nodes)))]
                node_sizes[this_begin_node] = parameters.pickup_node_scale
                node_sizes[this_end_node] = parameters.pickup_node_scale
                               
                networkx.draw(diGraphForPath, node_color = node_colors, node_size = node_sizes, width = 0.0, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = positions, with_labels = parameters.draw_graph_nodes_with_labels)

                #
                # Then draw the edges, setting a different color for each edge
                #

                colors = cm.rainbow(np.linspace(0, 1, len(this_shortest_path)))                
                edge_colors = [colors[i_step] for i_step in range(0, len(this_shortest_path) - 1)]
                networkx.draw_networkx_edges(diGraphForPath, pos = positions, edgelist = this_shortest_path_edges, edge_color = edge_colors, arrowsize = parameters.arrow_size)
                                
                ### END DEBUG

        
        #
        # Compute the mean harvestability over the given set of nodes
        #

        for n in list_of_nodes:

                this_node_attributes = the_graph.nodes[n]
        
                this_node_x_position = this_node_attributes['data'].x_position
                this_node_y_position = this_node_attributes['data'].y_position

                i_harvest = (np.round(-1.0*(this_node_y_position - harvestability_map_final_boundingbox_max_northing) / parameters.harvestability_map_pixel_size)).astype(int)
                j_harvest = (np.round((this_node_x_position - harvestability_map_final_boundingbox_min_easting) / parameters.harvestability_map_pixel_size)).astype(int)

                this_harvestability = harvestability_map[i_harvest, j_harvest]
                
                mean_harvestability = mean_harvestability + this_harvestability
                
        mean_harvestability = mean_harvestability / len(list_of_nodes)
                
        
        return the_distance_matrix, the_ground_damage_matrix, mean_harvestability



#
# Plot the harvest area real-life graph as decomposed into region
# graphs
#

def plot_decomposition(the_real_life_graph, pile_node_colors, pile_node_sizes, pile_positions_dictionary, decomposition, graph_axis_limits):

        number_of_region_graphs = len(decomposition)
        
        #
        # Ignore warnings for these plots
        #

        warnings.simplefilter("ignore")

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Decomposition of graph')


        #
        # First, plot all the region graphs into the same plot
        #
        
        #
        # Plot the nodes of the full real-life graph
        #

        networkx.draw(the_real_life_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)


        #
        # Plot the region graphs of the decomposition, each with a distinct color
        #
        
        i_region_graph = 0
        
        region_graph_colors = cm.rainbow(np.linspace(0, 1, number_of_region_graphs))
        
        for this_region_graph in decomposition:

                this_region_graph_graph = networkx.Graph()
                this_region_graph_graph.add_nodes_from(this_region_graph)
        
                networkx.draw(this_region_graph_graph, node_color = np.ones([len(list(this_region_graph_graph.nodes)), 4])*np.array(region_graph_colors[i_region_graph]), alpha = parameters.decomposition_alpha, node_size = parameters.decomposition_node_size, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, with_labels = parameters.draw_graph_nodes_with_labels)
        
                i_region_graph = i_region_graph + 1

        
        #
        # Add legend and set the axis limits to the same value as in
        # the plot of the full graph
        #

        decomposition_color_legend_handles = []

        i_region_graph = 1

        for region_graph_color in region_graph_colors:

                decomposition_color_legend_handles.append(mlines.Line2D([], [], color = region_graph_color, alpha = parameters.decomposition_alpha, marker = '.', linestyle = 'None', markersize = parameters.legend_region_graph_color_marker_size, label = 'Region graph ' + str(i_region_graph)))

                i_region_graph = i_region_graph + 1

        plt.legend(handles = decomposition_color_legend_handles, loc = 'upper right')
        plt.axis(graph_axis_limits)
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:
        
                plt.savefig('decomposition_of_real-life_graph.png')

        plt.close()
        

        #
        # Then, create a set of plots, one for each region graph
        #

        i_region_graph = 0

        for this_region_graph in decomposition:
        
                plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                plt.rc('font', size = parameters.default_font_size)
                plt.title('Region graph ' + str(i_region_graph + 1))

                #
                # First, plot the nodes of the full real-life graph
                #

                networkx.draw(the_real_life_graph, node_color = [c for c in 'k'*len(list(the_real_life_graph.nodes))], node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)


                #
                # Then, plot the nodes of the region graph in the designated color
                #
        
                this_region_graph_graph = networkx.Graph()
                this_region_graph_graph.add_nodes_from(this_region_graph)
        
                networkx.draw(this_region_graph_graph, node_color = np.ones([len(list(this_region_graph_graph.nodes)), 4])*np.array(region_graph_colors[i_region_graph]), alpha = parameters.decomposition_alpha, node_size = parameters.decomposition_node_size, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, with_labels = parameters.draw_graph_nodes_with_labels)
        
                #
                # Set axis limits to the same value as in the plot of
                # the full graph
                #

                plt.axis(graph_axis_limits)
                plt.gca().set_aspect('equal', adjustable = 'box')
        
                if parameters.save_figures_to_files:
        
                        plt.savefig('region_graph_' + str(i_region_graph + 1) + '.png')
        
                plt.close()
                        
                i_region_graph = i_region_graph + 1
                        
        
        #
        # Enable warnings again
        #

        warnings.resetwarnings()

        return


#
# Plot the best solution for a given subproblem on top of the
# real-life graph, first separately for each tour, and then all tours
# together
#

def plot_best_solution(the_real_life_graph, the_solution, i_subproblem, pile_node_colors, pile_node_sizes, pile_positions_dictionary, pile_color_legend_handles, graph_axis_limits, subproblem):

        #
        # Use a directed graph so that we can draw arrows
        #

        diGraphForPath = networkx.DiGraph()
        diGraphForPath.add_nodes_from(the_real_life_graph)

        
        #
        # Disable all warnings for this plot, to avoid unnecessary deprecation warnings
        #
        
        warnings.simplefilter("ignore")

        
        #
        # Plot each tour separately, if desired
        #
        
        if parameters.plot_each_tour_separately == True:

                for i_tour in range(1, the_solution.path.count(0)):

                        #
                        # Get the ith tour of the solution path
                        #
                        
                        this_tour_path = []
                        this_tour = 0
                        
                        for node in the_solution.path:

                                if node == 0:

                                        this_tour = this_tour + 1

                                if this_tour == i_tour:

                                        this_tour_path.append(node)

                        this_tour_path.append(0)

                        #
                        # Get the path edges for plotting
                        #
                        
                        this_tour_edges = [(this_tour_path[i], this_tour_path[i+1]) for i in range (0, len(this_tour_path) - 1)]

                        #
                        # Plot
                        #

                        # DEBUG
                        #fig, ax = plt.subplots(figsize = (parameters.figure_size_h, parameters.figure_size_v))

                        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
                        plt.rc('font', size = parameters.default_font_size)
                        plt.title('Tour number ' + str(i_tour) + ' best solution of subproblem ' + str(i_subproblem) + ' on real-life graph')

                        #
                        # Draw the real-life graph into the background first, along with the road connectivity. Color all nodes outside of this subproblem black.
                        #

                        pile_node_colors_for_this_subproblem = copy.deepcopy(pile_node_colors)

                        for node in range(0, len(pile_node_colors)):

                                if node > 0 and node not in subproblem:

                                        pile_node_colors_for_this_subproblem[node] = [0.0, 0.0, 0.0]
                        
                        # DEBUG
                        #networkx.draw(the_real_life_graph, node_color = pile_node_colors_for_this_subproblem, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels, ax = ax)

                        networkx.draw(the_real_life_graph, node_color = pile_node_colors_for_this_subproblem, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

                        #
                        # Then draw the path for this tour, setting the color of each edge to the color of the destination node of that edge
                        #

                        edge_colors = []
                        
                        for i_edge in range(0, len(this_tour_edges)):

                                next_node = this_tour_edges[i_edge][1]
                                edge_colors.append(pile_node_colors[next_node])
                
                        # DEBUG
                        #networkx.draw_networkx_edges(diGraphForPath, pos = pile_positions_dictionary, edgelist = this_tour_edges, edge_color = edge_colors, arrowsize = parameters.arrow_size, width = parameters.path_width, style = 'dashed', ax = ax)

                        networkx.draw_networkx_edges(diGraphForPath, pos = pile_positions_dictionary, edgelist = this_tour_edges, edge_color = edge_colors, arrowsize = parameters.arrow_size, width = parameters.path_width, style = 'dashed')

                        # DEBUG
                        #ax.tick_params(left = True, bottom = True, labelleft = True, labelbottom = True)

                        plt.legend(handles = pile_color_legend_handles, loc = 'lower left')
                        plt.axis(graph_axis_limits)
                        plt.gca().set_aspect('equal', adjustable = 'box')

                        # DEBUG
                        #plt.axis('on')
                        
                        if parameters.save_figures_to_files:
                        
                                plt.savefig('best_solution_tour_' + str(i_tour) + '_subproblem_' + str(i_subproblem) + '.png')
                        
                        plt.close()
                                

        #
        # In any case, plot the given solution with all tours into a single figure
        #

        #
        # Get the path edges for plotting
        #
        
        solution_path_edges = [(the_solution.path[i], the_solution.path[i+1]) for i in range (0, len(the_solution.path) - 1)]

        #
        # Plot
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.rc('font', size = parameters.default_font_size)
        plt.title('Best solution of subproblem ' + str(i_subproblem) + ' on real-life graph')

        #
        # Draw just the nodes first
        #                

        networkx.draw(the_real_life_graph, node_color = pile_node_colors, node_size = pile_node_sizes, font_color = 'k', font_size = parameters.graph_node_font_size, node_shape = 'o', pos = pile_positions_dictionary, edge_color = 'lightgrey', with_labels = parameters.draw_graph_nodes_with_labels)

        #
        # Then draw the edges, setting the color of each edge to the color of the destination node of that edge
        #
        
        edge_colors = []
        
        for i_edge in range(0, len(solution_path_edges)):
                
                next_node = solution_path_edges[i_edge][1]
                edge_colors.append(pile_node_colors[next_node])
        
        networkx.draw_networkx_edges(diGraphForPath, pos = pile_positions_dictionary, edgelist = solution_path_edges, edge_color = edge_colors, arrowsize = parameters.arrow_size, width = parameters.path_width)

        plt.legend(handles = pile_color_legend_handles, loc = 'upper right')
        plt.axis(graph_axis_limits)
        plt.gca().set_aspect('equal', adjustable = 'box')
        
        if parameters.save_figures_to_files:

                plt.savefig('best_solution_of_subproblem_' + str(i_subproblem) + '.png')

        plt.close()
        
        
        #
        # Enable warnings again
        #

        warnings.resetwarnings()

        return


#
# Plot total cost for best solution vs. optimization iteration
#

def plot_best_solution_cost_vs_iteration(i_subproblem, best_solution_cost_vs_iteration):

        #
        # Disable all warnings for this plot, to avoid unnecessary deprecation warnings
        #
        
        warnings.simplefilter("ignore")


        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))
        plt.title('Solution convergence in subproblem ' + str(i_subproblem))
        plt.plot(best_solution_cost_vs_iteration[:, 0], best_solution_cost_vs_iteration[:, 1], 'r-', linewidth=3)
        plt.xlabel('Iteration')
        plt.ylabel('Best solution cost')

        if parameters.save_figures_to_files:

                plt.savefig('best_solution_cost_vs_iteration_for_subproblem_' + str(i_subproblem) + '.png')

        plt.close()
        
        
        #
        # Enable warnings again
        #

        warnings.resetwarnings()

        return


#
# Plot the load composition vs. tour number for the best solution so
# far
#

def plot_load_composition_vs_tour(best_solution_so_far, wood_matrix, wood_type_dictionary, new_wood_type_dictionary, new_wood_type_color_dictionary, i_subproblem):

        #
        # Disable all warnings for this plot, to avoid unnecessary deprecation warnings
        #
        
        warnings.simplefilter("ignore")

        
        #
        # Separate the total path into individual tours
        #
        
        individual_tours = []
        this_tour = []
        
        for node in best_solution_so_far.path:
        
                if node != 0:

                        this_tour.append(node)
                
                else:

                        individual_tours.append(this_tour)
                        this_tour = []

        #
        # Discard the first tour, as this is equal to just []
        #
        
        individual_tours = individual_tours[1:]
        number_of_tours = len(individual_tours)
        
        wood_types = np.unique(wood_matrix[1:, 0])
        amounts_of_each_wood_type = np.zeros([len(wood_type_dictionary), number_of_tours])

        for i_tour in range(0, number_of_tours):

                for node in individual_tours[i_tour]:

                        wood_type = wood_matrix[node, 0]                        

                        #
                        # The matrix indeces correspond to the (new) wood types as "index 0 : wood type 1", "index 1 : wood type 2", etc.
                        #
                        
                        amounts_of_each_wood_type[int(wood_type) - 1, i_tour] = amounts_of_each_wood_type[int(wood_type) - 1, i_tour] + wood_matrix[node, 1]
                

        #
        # Plot the data as a stacked bar diagram
        #

        plt.figure(figsize = (parameters.figure_size_h, parameters.figure_size_v))

        tour_numbers = np.arange(1, number_of_tours + 1)
        plot_handles = []
        legend_strings = []

        for wt in range(0, len(wood_type_dictionary)):

                if wt == 0:
                        
                        this_bottom = np.zeros(number_of_tours)
                else:
                        
                        this_bottom = np.sum(amounts_of_each_wood_type[0:wt, :], axis = 0)
                
                h = plt.bar(tour_numbers, amounts_of_each_wood_type[wt, :], width = parameters.composition_bar_chart_width, bottom = this_bottom, color = new_wood_type_color_dictionary[wt + 1])
        
                plot_handles.append(h)
                legend_strings.append('Assortment ' + str(new_wood_type_dictionary[wt + 1]))

        plt.legend(plot_handles, legend_strings)
        plt.title('Composition of load in subproblem ' + str(i_subproblem))
        plt.xlabel('Tour number')
        plt.xticks(np.arange(0, number_of_tours) + 1)
        plt.ylabel('Volume (m$^3$)')        
        
        if parameters.save_figures_to_files:

                plt.savefig('best_solution_load_composition_in_subproblem_' + str(i_subproblem) + '.png')

        plt.close()
        

        #
        # Enable warnings again
        #

        warnings.resetwarnings()

        return

