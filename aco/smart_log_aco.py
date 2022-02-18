#!/usr/bin/python3
#
# -*- coding: utf-8 -*-
#
# smart_log_aco.py
#
# Heuristically solves a multi-objective vehicle routing problem for
# optimizing forest forwarder work in a cut-to-length harvesting
# operation.
#
# Copyright 2022 Natural Resources Institute Finland
#
# This program is distributed under the terms of the GNU Lesser
# General Public License version 3.0 or any later version.
#


import sys
import os

#
# Try to import parameters.py from the local directory first. To accomplish this,
# insert the current working directory to the head of the system path list.
#

sys.path.insert(0, os.getcwd())

import parameters

#
# Revert the path back to look for other modules elsewhere than in the
# current working directory.
#

sys.path = sys.path[1:]

#
# Then import the rest of the modules.
#

import smart_log_classes_aco
import functions
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import time


#
# This is required for CSC batch jobs
#

matplotlib.use("agg")


#
# Usage
#

if len(sys.argv) < 3:

        print("Usage: smart_log_aco.py [wood pile input data file] [forest machine position input data file]")
        exit(0)


#
# Save out the parameters.py file that was used here
#

parameters_full_path = os.path.abspath(parameters.__file__)

print("")
print("The module parameters.py was loaded from the file %s. Now copying this file to parameters.out..." % parameters_full_path)
os.system("cp " + parameters_full_path + " parameters.out")
print("Done.")


#
# Measure the total wall-clock time consumed by the program
#

t_begin_program = time.time()


#
# Assign input parameters
#

wood_pile_input_data_file = sys.argv[1]
forwarder_position_input_data_file = sys.argv[2]


#
# Print out the value of each parameter registered in parameters.py
#

print("")
print("Found the following parameter values in parameters.py:")
print("")

print("random_number_generator_seed =", parameters.random_number_generator_seed)
print("harvestability_map_file =", parameters.harvestability_map_file)
print("harvestability_map_ulc_pixel_center_x_coordinate =", parameters.harvestability_map_ulc_pixel_center_x_coordinate)
print("harvestability_map_ulc_pixel_center_y_coordinate =", parameters.harvestability_map_ulc_pixel_center_y_coordinate)
print("harvestability_map_pixel_size =", parameters.harvestability_map_pixel_size) 
print("x_roadside_position =", parameters.x_roadside_position)
print("y_roadside_position =", parameters.y_roadside_position)
print("path_to_visualize =", parameters.path_to_visualize)
print("path_to_visualize_axis_limits =", parameters.path_to_visualize_axis_limits)
print("forwarder_initial_capacity =", parameters.forwarder_initial_capacity)
print("load_full_tolerance =", parameters.load_full_tolerance)
print("forwarder_speed =", parameters.forwarder_speed)
print("save_figures_to_files =", parameters.save_figures_to_files)
print("plot_each_tour_separately =", parameters.plot_each_tour_separately)
print("save_joining_and_merging_figures_to_files =", parameters.save_joining_and_merging_figures_to_files)
print("save_individual_subgraph_plots_to_files =", parameters.save_individual_subgraph_plots_to_files)
print("pickup_node_marker_scaling_factor =", parameters.pickup_node_marker_scaling_factor)
print("default_pile_marker_size =", parameters.default_pile_marker_size)
print("legend_wood_pile_marker_size =", parameters.legend_wood_pile_marker_size)
print("roadside_node_scale =", parameters.roadside_node_scale)
print("road_point_markersize =", parameters.road_point_markersize)
print("road_point_marker_edge_width =", parameters.road_point_marker_edge_width)
print("raw_forwarder_data_markersize =", parameters.raw_forwarder_data_markersize)
print("raw_forwarder_position_point_edge_width =", parameters.raw_forwarder_position_point_edge_width)
print("arrow_size =", parameters.arrow_size)
print("path_width =", parameters.path_width)
print("road_strip_width =", parameters.road_strip_width)
print("node_size_in_road_strip_graph =", parameters.node_size_in_road_strip_graph)
print("fragment_path_width =", parameters.fragment_path_width)
print("fragment_graph_node_size =", parameters.fragment_graph_node_size)
print("subgraph_path_width =", parameters.subgraph_path_width)
print("subgraph_cm_marker_size =", parameters.subgraph_cm_marker_size)
print("subgraph_cm_marker_edge_width =", parameters.subgraph_cm_marker_edge_width)
print("node_size_in_before_and_after_joining_plots =", parameters.node_size_in_before_and_after_joining_plots)
print("decomposition_alpha =", parameters.decomposition_alpha)
print("legend_region_graph_color_marker_size =", parameters.legend_region_graph_color_marker_size)
print("figure_size_h =", parameters.figure_size_h)
print("figure_size_v =", parameters.figure_size_v)
print("margin_for_pile_and_road_plots =", parameters.margin_for_pile_and_road_plots)
print("n_bins_wood_type_histogram =", parameters.n_bins_wood_type_histogram)
print("pile_volume_histogram_bin_step =", parameters.pile_volume_histogram_bin_step) 
print("pile_number_histogram_bin_step =", parameters.pile_number_histogram_bin_step)
print("pile_extent_histogram_bin_step =", parameters.pile_extent_histogram_bin_step)
print("composition_bar_chart_width =", parameters.composition_bar_chart_width) 
print("small_font_size =", parameters.small_font_size) 
print("default_font_size =", parameters.default_font_size) 
print("graph_node_font_size =", parameters.graph_node_font_size)
print("draw_graph_nodes_with_labels =", parameters.draw_graph_nodes_with_labels)
print("title_pad =", parameters.title_pad)
print("maximum_number_of_nodes_per_region_graph =", parameters.maximum_number_of_nodes_per_region_graph)
print("minimum_number_of_nodes_per_subgraph =", parameters.minimum_number_of_nodes_per_subgraph)
print("decomposition_node_size =", parameters.decomposition_node_size)
print("road_point_identification_window_size =", parameters.road_point_identification_window_size)
print("road_point_identification_min_n_points_per_window =", parameters.road_point_identification_min_n_points_per_window)
print("n_road_points_for_pile_placement =", parameters.n_road_points_for_pile_placement)
print("road_spline_parameter_step =", parameters.road_spline_parameter_step)
print("pile_connectivity_initial_sector_range =", parameters.pile_connectivity_initial_sector_range)
print("pile_connectivity_sector_range_increase =", parameters.pile_connectivity_sector_range_increase)
print("pile_connectivity_max_distance_between_neighbors =", parameters.pile_connectivity_max_distance_between_neighbors)
print("pile_connectivity_max_sector =", parameters.pile_connectivity_max_sector)
print("t_load_a =", parameters.t_load_a)
print("t_load_b =", parameters.t_load_b)
print("t_load_c =", parameters.t_load_c)
print("t_unload_a =", parameters.t_unload_a)
print("t_unload_b =", parameters.t_unload_b)
print("t_unload_c =", parameters.t_unload_c)
print("minimum_cost_in_distance_matrix =", parameters.minimum_cost_in_distance_matrix)
print("ground_damage_model_v0 =", parameters.ground_damage_model_v0)
print("cap_on_effective_number_of_unique_wood_types =", parameters.cap_on_effective_number_of_unique_wood_types)
print("w_length =", parameters.w_length)
print("w_duration =", parameters.w_duration)
print("w_ground_damage =", parameters.w_ground_damage)
print("tau_0 =", parameters.tau_0)
print("rho =", parameters.rho)
print("alpha =", parameters.alpha)
print("beta =", parameters.beta)
print("number_of_ants =", parameters.number_of_ants)
print("number_of_candidates =", parameters.number_of_candidates)
print("max_iterations =", parameters.max_iterations)


#
# Seed both the numpy and non-numpy random number generators. Numpy
# used for all cases where randomness in numbers is needed. The
# non-numpy generator is used only for shuffling colors in the plots.
#

print("")

if parameters.random_number_generator_seed >= 0:

        this_seed = parameters.random_number_generator_seed

        np.random.seed(seed = this_seed)
        random.seed(a = this_seed)
        
        print("Seeding the random number generators with seed %d" % this_seed)
        
else:

        np.random.seed(seed = None)
        random.seed(a = None)
        
        print("Seeding the random number generators from system data.")

        
#
# Create a real-world graph of the harvest area
#

t_begin = time.time()

the_real_life_graph, pile_node_sizes, pile_node_colors, pile_color_legend_handles, wood_type_color_dictionary, wood_type_dictionary, graph_axis_limits, pile_positions_dictionary = functions.construct_real_life_graph_from_input_files(wood_pile_input_data_file, forwarder_position_input_data_file)

t_end = time.time()
time_for_creating_real_life_graph = t_end - t_begin


#
# If a path was given for visualization in parameters.py, visualize
# the path and then exit.
#

if parameters.path_to_visualize != []:

        print("")
        print("Now visualizing the path given in parameters.py. NB! The set of nodes of the path will be treated as the subproblem in the plot.")
        
        solution_to_visualize = smart_log_classes_aco.solution(parameters.path_to_visualize, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf)
        functions.plot_best_solution(the_real_life_graph, solution_to_visualize, 1, pile_node_colors, pile_node_sizes, pile_positions_dictionary, pile_color_legend_handles, parameters.path_to_visualize_axis_limits, set(parameters.path_to_visualize))

        plt.show()
        
        print("Exiting.")
        print("")

        exit(0)



#
# Create a new wood type color dictionary of the format
#
# {new wood type : color}
#
# as in this module we deal with the new labels
#
# (1, 2, 3, ...)
#
# instead of whatever was given in the original data and which then
# was used in the module functions.py.
#

new_wood_type_color_dictionary = {}

for wood_type in wood_type_dictionary:
        
        new_wood_type_color_dictionary[wood_type_dictionary[wood_type]] = wood_type_color_dictionary[wood_type]


#
# Also, create a new wood type dictionary of the format
#
# {new label : original label}
#
# as we deal with the new labels
#
# (1, 2, 3, ...) in this module.
#

new_wood_type_dictionary = {}

for wood_type in wood_type_dictionary:
        
        new_wood_type_dictionary[wood_type_dictionary[wood_type]] = wood_type


#
# Get a static harvestability map encompassing the entire harvest
# area. This will be used for computing the cost matrix for each
# individual subproblem later on.
#

#
# Find the min-max values of northing and easting for the piles in the harvest area
#

print("")
print("Now finding the minimum and maximum values of pile positions in the real-life graph...")

harvest_area_boundingbox_min_easting = np.Inf
harvest_area_boundingbox_max_easting = 0.0
harvest_area_boundingbox_min_northing = np.Inf
harvest_area_boundingbox_max_northing = 0.0

for node in pile_positions_dictionary:
        
        this_pile_easting = pile_positions_dictionary[node][0]
        this_pile_northing = pile_positions_dictionary[node][1]

        if this_pile_easting < harvest_area_boundingbox_min_easting:
                
                harvest_area_boundingbox_min_easting = this_pile_easting
        
        if this_pile_easting > harvest_area_boundingbox_max_easting:
                
                harvest_area_boundingbox_max_easting = this_pile_easting

        if this_pile_northing < harvest_area_boundingbox_min_northing:
                
                harvest_area_boundingbox_min_northing = this_pile_northing
        
        if this_pile_northing > harvest_area_boundingbox_max_northing:
                
                harvest_area_boundingbox_max_northing = this_pile_northing


print("Done. Node positions run from easting %f m to %f m and northing %f m to %f m" % (harvest_area_boundingbox_min_easting, harvest_area_boundingbox_max_easting, harvest_area_boundingbox_min_northing, harvest_area_boundingbox_max_northing))


#
# Then, get the harvestability map
#

print("")
print("Now getting the harvestability map for the entire harvest area.")

harvestability_map, harvestability_map_final_boundingbox_min_easting, harvestability_map_final_boundingbox_max_easting, harvestability_map_final_boundingbox_min_northing, harvestability_map_final_boundingbox_max_northing = functions.get_harvestability_map(the_real_life_graph, pile_node_sizes, pile_positions_dictionary, graph_axis_limits, harvest_area_boundingbox_min_easting, harvest_area_boundingbox_max_easting, harvest_area_boundingbox_min_northing, harvest_area_boundingbox_max_northing)

print("Done.")


#
# Split the full graph into smaller graphs. The optimization
# problem will be solved separately for each of these.
#

print("")
print("Now starting decomposition procedure for the real-life graph.")

t_begin = time.time()
        
decomposition = functions.decompose_graph(the_real_life_graph, pile_node_colors, pile_node_sizes, pile_color_legend_handles, pile_positions_dictionary, graph_axis_limits)

t_end = time.time()
time_for_creating_decomposition = t_end - t_begin


#
# Get some statistics for the decomposition and print these out
#

number_of_region_graphs = len(decomposition)

sizes_of_region_graphs = []

for this_region_graph in decomposition:

        sizes_of_region_graphs.append(len(this_region_graph))

sizes_of_region_graphs = np.array(sizes_of_region_graphs)

minimum_size_of_region_graph = np.min(sizes_of_region_graphs)
maximum_size_of_region_graph = np.max(sizes_of_region_graphs)
mean_size_of_region_graph = np.mean(sizes_of_region_graphs)
std_of_region_graph_size = np.std(sizes_of_region_graphs)

total_number_of_nodes = len(list(the_real_life_graph.nodes))
total_number_of_pickup_nodes = total_number_of_nodes - 1

print("")
print("Decomposition procedure complete.")
print("The decomposition consists of a total of %d region graphs with size ranging from %d to %d nodes, with a mean of %f nodes and a standard deviation of %f nodes." % (number_of_region_graphs, minimum_size_of_region_graph, maximum_size_of_region_graph, mean_size_of_region_graph, std_of_region_graph_size))
print("")
print("The total number of nodes in the entire system is %d, of which %d are pickup nodes." % (total_number_of_nodes, total_number_of_pickup_nodes))


#
# Plot the decomposition
#

functions.plot_decomposition(the_real_life_graph, pile_node_colors, pile_node_sizes, pile_positions_dictionary, decomposition, graph_axis_limits)


#
# Loop over the region graphs of the decomposition, solving the optimization problem separately
# for each.
#

print("")
print("*** Now starting loop over subproblems ***")

i_subproblem = 1

total_time_for_optimization = 0.0
total_time_for_computing_cost_matrices = 0.0

total_weighted_cost = 0.0
total_length = 0.0
total_duration = 0.0
total_ground_damage = 0.0

total_number_of_tours = 0
total_number_of_piles = 0
total_volume_of_wood_forwarded = 0.0


for subproblem in decomposition:

        
        #
        # Print out essential parameter values
        #
        
        print("")
        print("Using the following parameter values as read from parameters.py:")
        print("")
        print("tau_0 = %f" % parameters.tau_0)
        print("rho = %f" % parameters.rho)
        print("alpha = %f" % parameters.alpha)
        print("beta = %f" % parameters.beta)
        print("number_of_ants = %d" % parameters.number_of_ants)
        print("number_of_candidates = %d" % parameters.number_of_candidates)
        print("max_iterations = %d" % parameters.max_iterations)
        print("forwarder_initial_capacity = %f" % parameters.forwarder_initial_capacity)
        print("forwarder_speed = %f" % parameters.forwarder_speed)
        print("w_length = %f" % parameters.w_length)
        print("w_duration = %f" % parameters.w_duration)
        print("w_ground_damage = %f" % parameters.w_ground_damage)
        
        
        #
        # Keep track of the mean values of the components of the
        # heuristic desirability
        #
        
        sum_of_distance_cost_ij = 0.0
        sum_of_duration_cost_ij = 0.0
        sum_of_ground_damage_cost_ij = 0.0

        n_distance_cost_ij = 0
        n_duration_cost_ij = 0
        n_ground_damage_cost_ij = 0
        
        
        print("")
        print("===> Now starting to process subproblem number %d <===" % i_subproblem)

        t_begin_subproblem = time.time()
        
        #
        # Add the roadside node to the subproblem
        #

        subproblem.insert(0, 0)

        #
        # Store the total number of nodes and pickup nodes for future
        # use
        #
        
        number_of_nodes_in_subproblem = len(subproblem)
        number_of_pickup_nodes_in_subproblem = number_of_nodes_in_subproblem - 1

        print("")
        print("This subproblem consists of %d pickup nodes and one roadside node." % number_of_pickup_nodes_in_subproblem)

        
        #
        # If desired, set the number of ants equal to the number of pickup
        # nodes. Otherwise, cap the number of ants to the number of pickup
        # nodes.
        #

        if parameters.number_of_ants < 0:        

                print("")
                print("Setting number of ants to number of nodes in subproblem")
                this_subproblem_number_of_ants = number_of_pickup_nodes_in_subproblem

        else:

                this_subproblem_number_of_ants = np.min([number_of_pickup_nodes_in_subproblem, parameters.number_of_ants])


        print("")
        print("The number of ants used in this subproblem is %d" % this_subproblem_number_of_ants)
                
        #
        # Based on shortest paths in the real-world graph, compute the
        # (generally sparse) distance matrix D_ij for this
        # subproblem. Also, compute the ground damage matrix G_ij for
        # this problem.
        #
        
        t_begin = time.time()

        print("")
        print("Computing the static cost matrices from scratch...")
        
        distance_matrix, ground_damage_matrix, mean_harvestability = functions.compute_static_cost_matrices(the_real_life_graph, subproblem, harvestability_map, harvestability_map_final_boundingbox_min_easting, harvestability_map_final_boundingbox_max_easting, harvestability_map_final_boundingbox_min_northing, harvestability_map_final_boundingbox_max_northing)

        print("Done.")

        #
        # Compute statistics on distance matrix
        #

        masked_distance_matrix = np.ma.masked_where(distance_matrix == np.inf, distance_matrix) 

        min_distance_matrix = np.min(masked_distance_matrix)
        mean_distance_matrix = np.mean(masked_distance_matrix)
        max_distance_matrix = np.max(masked_distance_matrix)
        std_distance_matrix = np.std(masked_distance_matrix)
        
        print("")
        print("Distance matrix values go from %f m to %f m with a mean of %f m and a std of %f m" % (min_distance_matrix, max_distance_matrix, mean_distance_matrix, std_distance_matrix))
        
        print("")
        print("Found a mean harvestability of %f over the nodes in this subproblem." % mean_harvestability)
        
        
        t_end = time.time()
        time_for_computing_cost_matrices_for_this_subproblem = t_end - t_begin
        total_time_for_computing_cost_matrices = total_time_for_computing_cost_matrices + time_for_computing_cost_matrices_for_this_subproblem
        
        
        #
        # Create the wood pile data matrix. The matrix holds the wood
        # type (assortment) and wood amount (m**3) for each pickup
        # node. Matrix index is equal to node name.
        #
        
        wood_matrix = np.zeros([total_number_of_nodes, 2])
        
        for node in subproblem:

                #
                # Only consider nodes other than the roadside node
                #
                
                if node != 0:
                
                        this_wood_type = the_real_life_graph.nodes[node]['data'].wood_pile.type_of_wood
                        this_wood_amount = the_real_life_graph.nodes[node]['data'].wood_pile.amount_of_wood
                
                        wood_matrix[node, 0] = this_wood_type
                        wood_matrix[node, 1] = this_wood_amount


        #
        # Create the candidate list for each pickup node. The
        # candidates for a given node are the number_of_candidates
        # closest nodes by distance D_ij.
        #
        
        print("")
        print("Now forming the candidate lists...")
        
        #
        # Set number of candidates per node to the default value, if
        # desired. Otherwise, cap the number of candidates to the
        # number of pickup nodes.
        #

        if parameters.number_of_candidates < 0:

                print("Setting number of candidates to number of nodes in subproblem / 4")
                this_subproblem_number_of_candidates = int(number_of_nodes_in_subproblem / 4.0)

        else:
                
                this_subproblem_number_of_candidates = np.min([number_of_pickup_nodes_in_subproblem, parameters.number_of_candidates])
                

        print("")
        print("The number of candidates used in this subproblem is %d " % this_subproblem_number_of_candidates)
                
        candidate_lists = {}

        #
        # Construct candidate list for each pickup node
        #
        
        for node in subproblem[1:]:

                #
                # Find the closest nodes to node i_node
                #
        
                the_closest_nodes = list(np.argsort(distance_matrix[node, :])[0 : this_subproblem_number_of_candidates])
                
                #
                # Only nodes within this subproblem are allowed in the
                # candidate list. Enforce this here.
                #
                
                for candidate_node in the_closest_nodes:

                        if candidate_node not in subproblem:

                                print("Warning! Node %d was found in candidate list of node %d but is not included in this subproblem. Now removing it from the candidate list." % (candidate_node, node))
                                the_closest_nodes.remove(candidate_node)

                #
                # Save the candidate list for this node
                #
                                
                candidate_lists[node] = the_closest_nodes

                
        #
        # Add all pickup nodes to the candidate list of the roadside
        # node
        #
        
        candidate_lists[0] = [node for node in subproblem[1:]]

        print("")
        print("Done.")

        
        #
        # Compute the total and mean volume of piles for this
        # subproblem
        #
        
        total_volume_of_piles = np.sum(wood_matrix[:, 1])
        mean_volume_of_piles = total_volume_of_piles / number_of_pickup_nodes_in_subproblem
        
        print("")
        print("Total volume of wood piles is %f m**3, mean volume is %f m**3" % (total_volume_of_piles, mean_volume_of_piles))

        
        #
        # Also, get the number of unique assortments for this
        # subproblem
        #
        
        number_of_unique_wood_types_for_this_subproblem = len(np.unique(wood_matrix[subproblem[1:], 0]))

        print("The number of unique assortments is %d" % number_of_unique_wood_types_for_this_subproblem)
        
        #
        # Use the distance matrix, the ground damage matrix, and
        # forwarder speed along with loading and unloading times to
        # compute an estimate for the length, duration, and ground
        # damage of a typical path, i.e., a solution for this
        # subproblem. These are needed for standardizing the objective
        # components in the total cost of the path.
        #
        # In addition, compute the mean of the ground damage matrix,
        # the mean time of traversing an arc, and the mean damage
        # created upon traversing an arc. These are needed for
        # standardizing the heuristic desirability in the computation
        # of p_ij.
        #

        print("")
        print("Now estimating typical values for total distance, total duration, and total ground damage of a path in this subproblem.")
        print("These will be used for standardizing the objective components in the cost function approximately to the neighborhood of 1.0")
        print("")
        print("Also, estimating typical values for the components of the heuristic desirability in p_ij.")
        print("These will be used for standardizing the components approximately to the neighborhood of 1.0")

        
        #
        # Compute statistics on the mean nearest-neighbor distance of loading nodes
        #

        nn_distances = []
        
        for node in subproblem[1:]:
 
                this_nn_distance = distance_matrix[node, candidate_lists[node][0]]
                nn_distances.append(this_nn_distance)

        nn_distances = np.array(nn_distances)

        min_nn_distance = np.min(nn_distances)
        mean_nn_distance = np.mean(nn_distances)
        max_nn_distance = np.max(nn_distances)
        std_of_nn_distance = np.std(nn_distances)
        
        print("")
        print("Nearest-neighbor distance of loading nodes ranges from %f m to %f m with a mean of %f m and std of %f m" % (min_nn_distance, max_nn_distance, mean_nn_distance, std_of_nn_distance))
        
        #
        # Then, compute the mean distance from the unloading node to the loading nodes
        #

        x_mean_loading_nodes = 0.0
        y_mean_loading_nodes = 0.0
        
        for node in subproblem[1:]:
        
                this_node_attributes = the_real_life_graph.nodes[node]
        
                this_node_x_position = this_node_attributes['data'].x_position
                this_node_y_position = this_node_attributes['data'].y_position

                x_mean_loading_nodes = x_mean_loading_nodes + this_node_x_position
                y_mean_loading_nodes = y_mean_loading_nodes + this_node_y_position

        x_mean_loading_nodes = x_mean_loading_nodes / number_of_pickup_nodes_in_subproblem
        y_mean_loading_nodes = y_mean_loading_nodes / number_of_pickup_nodes_in_subproblem        

        mean_distance_from_unloading_node_to_loading_nodes = np.sqrt((parameters.x_roadside_position - x_mean_loading_nodes)**2 + (parameters.y_roadside_position - y_mean_loading_nodes)**2)

        #
        # From these, estimate the mean distance of traversing an arc
        # in a typical path
        #
        
        mean_of_traversal_distance_between_nodes = 2.0*mean_nn_distance
        total_traversal_distance_between_unloading_node_and_loading_nodes = 2.0 * 2.0 * mean_distance_from_unloading_node_to_loading_nodes * total_volume_of_piles / parameters.forwarder_initial_capacity
        
        #
        # From these, compute an estimate for the typical path length
        #

        estimated_mean_length_of_path = mean_of_traversal_distance_between_nodes*number_of_nodes_in_subproblem + total_traversal_distance_between_unloading_node_and_loading_nodes

        #
        # Furthermore, compute an estimate for the distance cost in the heuristic desirability
        #

        estimate_of_distance_cost_in_heuristic_desirability = 0.5*np.mean(distance_matrix[np.where(distance_matrix != np.Inf)])


        
        #
        # Estimate the typical total duration of a path
        #

        effective_number_of_unique_wood_types_for_this_subproblem = np.min([parameters.cap_on_effective_number_of_unique_wood_types, number_of_unique_wood_types_for_this_subproblem])

        print("Effective number of unique wood types for this subproblem is %d" % effective_number_of_unique_wood_types_for_this_subproblem)
        
        time_of_loading_one_m3 = parameters.t_load_a * np.power(effective_number_of_unique_wood_types_for_this_subproblem, 2) + parameters.t_load_b * effective_number_of_unique_wood_types_for_this_subproblem + parameters.t_load_c
        time_of_loading_one_m3 = time_of_loading_one_m3*60.0
        
        time_of_unloading_one_m3 = parameters.t_unload_a * np.power(effective_number_of_unique_wood_types_for_this_subproblem, 2) + parameters.t_unload_b * effective_number_of_unique_wood_types_for_this_subproblem + parameters.t_unload_c
        time_of_unloading_one_m3 = time_of_unloading_one_m3*60.0

        estimated_mean_time_of_loading_plus_unloading = (time_of_loading_one_m3 + time_of_unloading_one_m3)*total_volume_of_piles
        
        estimated_mean_duration_of_path = estimated_mean_length_of_path / parameters.forwarder_speed + estimated_mean_time_of_loading_plus_unloading
        
        #
        # Furthermore, compute an estimate for the duration cost in the heuristic desirability
        #

        estimate_of_duration_cost_in_heuristic_desirability = estimate_of_distance_cost_in_heuristic_desirability / parameters.forwarder_speed


                
        #
        # Estimate the mean ground damage created when traversing an
        # arc in a typical path
        #
        
        mean_of_traversal_ground_damage_between_nodes = mean_harvestability*mean_of_traversal_distance_between_nodes*(1.0 + 0.5*parameters.forwarder_initial_capacity / parameters.ground_damage_model_v0)
        total_traversal_ground_damage_between_unloading_node_and_loading_nodes = mean_harvestability*total_traversal_distance_between_unloading_node_and_loading_nodes*(1.0 + 0.5*parameters.forwarder_initial_capacity / parameters.ground_damage_model_v0)
        
        #
        # Estimate the typical ground damage of a given path
        #

        estimated_mean_ground_damage_of_path = mean_of_traversal_ground_damage_between_nodes*number_of_nodes_in_subproblem + total_traversal_ground_damage_between_unloading_node_and_loading_nodes
        
        #
        # Furthermore, compute an estimate for the ground damage cost in the heuristic desirability
        #

        estimate_of_ground_damage_cost_in_heuristic_desirability = estimate_of_distance_cost_in_heuristic_desirability*mean_harvestability*(1.0 + 0.5*parameters.forwarder_initial_capacity / parameters.ground_damage_model_v0)


        
        #
        # From the above, compute an estimate for the total cost of a
        # typical path in this subproblem
        #
        
        estimated_total_weighted_cost_of_path = parameters.w_length*1.0 + parameters.w_duration*1.0 + parameters.w_ground_damage*1.0


        
        print("")
        print("Done. For this subproblem, estimates of the mean values of objectives are as follows:")
        print("")
        print("Mean length of path: %f" % estimated_mean_length_of_path)
        print("Mean duration of path: %f" % estimated_mean_duration_of_path)
        print("Mean ground damage of path: %f" % estimated_mean_ground_damage_of_path)
        print("")
        print("Estimates of the mean values of the components in the heuristic desirability are as follows:")
        print("")
        print("Distance cost: %f" % estimate_of_distance_cost_in_heuristic_desirability)
        print("Duration cost: %f" % estimate_of_duration_cost_in_heuristic_desirability)
        print("Ground damage cost: %f" % estimate_of_ground_damage_cost_in_heuristic_desirability)
        print("")
        print("An estimate for the typical weighted cost of a path in this subproblem is %f*1.0 + %f*1.0 + %f*1.0 = %f" % (parameters.w_length, parameters.w_duration, parameters.w_ground_damage, estimated_total_weighted_cost_of_path))
        
        print("")
        print("Done.")

        
        #
        # Create the pheromone matrix, where the matrix element ij
        # gives the pheromone concentration on the path between nodes
        # i and j. If desired, initialize the matrix elements to
        #
        # tau_0 = 1 / (number_of_nodes * estimated_total_weighted_cost_of_path)
        #
        # Otherwise, use the value of tau_0 given in parameters.py.
        #

        print("")
        print("Now creating the pheromone matrix...")
        
        if parameters.tau_0 < 0.0:

                print("Setting tau_0 to 1 / (number of nodes in subproblem * estimated total weighted cost of path)")
                this_subproblem_tau_0 = 1.0 / (number_of_nodes_in_subproblem * estimated_total_weighted_cost_of_path)

        else:

                this_subproblem_tau_0 = parameters.tau_0

        print("")
        print("The value of tau_0 used in this subproblem is %f" % this_subproblem_tau_0)
        
        pheromone_matrix = np.ones([total_number_of_nodes, total_number_of_nodes])*this_subproblem_tau_0

        print("")
        print("Done.")

        
        #
        # Optimize the path using ACO. The goal is to find the
        # minimum cost path for the forwarder so that:
        #
        # - Each pickup node is visited exactly once
        #
        # - Each individual tour begins and ends at the
        #   roadside node
        #
        # - For each individual tour, wood piles are picked up
        #   until the forwarder capacity is exceeded
        #
        # Pheromone updating is done as follows (following the
        # original ACS):
        #
        # 1. When any given ant is constructing its tour,
        # update the edge ij that a given ant selects by
        #
        #   tau_ij <- (1 - rho) * tau_ij + rho * tau_0
        #
        #   where tau_0 is the initial value of pheromone
        #   concentration on each edge.
        #
        # 2. Every time the whole set of ants has completed a
        # tour, update the edges belonging to the best tour so
        # far (with cost C) according to
        #
        #   tau_ij <- (1 - rho) * tau_ij + rho / C
        #
        #   where rho is a parameter governing pheromone decay
        #
        # Cost is a weighted sum of distance traveled, time
        # used, and damage inflicted on the ground.
        #

        t_begin_optimization_of_this_subproblem = time.time()

        time_for_computing_solution_cost_for_this_subproblem = 0.0

        best_solution_so_far = smart_log_classes_aco.solution([], np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf)

        best_solution_cost_vs_iteration = np.zeros([parameters.max_iterations, 2])


        #
        # Start the iteration
        #

        print("")
        print("Starting iteration.")
        print("")

        for i_iteration in range(0, parameters.max_iterations):

                for i_ant in range(0, this_subproblem_number_of_ants):

                        #
                        # Use this ant to pick up all the wood
                        # piles, i.e., find a full set of
                        # tours for visiting all the loading
                        # nodes, i.e., construct a solution
                        #                

                        unvisited_pickup_nodes = [node for node in subproblem[1:]]

                        this_path = [0]
                        current_node = 0

                        #
                        # Keep track of the number of the tour
                        # that the forwarder is currently
                        # performing
                        #

                        i_tour = 0

                        while len(unvisited_pickup_nodes) > 0:


                                #
                                # Keep track of the volume of the current forwarder load
                                #

                                forwarder_load_volume = 0.0
                                i_tour = i_tour + 1

                                #
                                # Construct the next tour in this solution
                                #

                                while parameters.forwarder_initial_capacity - forwarder_load_volume >= parameters.load_full_tolerance and len(unvisited_pickup_nodes) > 0:

                                        #
                                        # For the first transition from the unloading node, set the ant on a randomly
                                        # selected node. After that, choose the node using the candidate list until it
                                        # is depleted.
                                        #

                                        if len(this_path) == 1:

                                                next_node = np.random.choice(subproblem[1:])

                                        else:

                                                unvisited_candidates = list((set(candidate_lists[current_node])).intersection(set(unvisited_pickup_nodes)))

                                                if len(unvisited_candidates) > 0:

                                                        distance_cost_ij = distance_matrix[current_node, unvisited_candidates] / estimate_of_distance_cost_in_heuristic_desirability
                                                        duration_cost_ij = distance_matrix[current_node, unvisited_candidates] / parameters.forwarder_speed / estimate_of_duration_cost_in_heuristic_desirability
                                                        ground_damage_cost_ij = ground_damage_matrix[current_node, unvisited_candidates]*(1.0 + forwarder_load_volume / parameters.ground_damage_model_v0) / estimate_of_ground_damage_cost_in_heuristic_desirability

                                                        total_cost_ij = parameters.w_length*distance_cost_ij + parameters.w_duration*duration_cost_ij + parameters.w_ground_damage*ground_damage_cost_ij

                                                        tau_ij = np.power(pheromone_matrix[current_node, unvisited_candidates], parameters.alpha)
                                                        eta_ij = np.power(1.0 / total_cost_ij, parameters.beta)

                                                        p_ij = tau_ij*eta_ij / np.sum(tau_ij*eta_ij)

                                                        i_next_node = np.random.choice(np.arange(0, len(unvisited_candidates)), p = p_ij)
                                                        next_node = unvisited_candidates[i_next_node]


                                                        #
                                                        # Update sums for the components of desirability
                                                        #

                                                        sum_of_distance_cost_ij = sum_of_distance_cost_ij + np.sum(distance_cost_ij)
                                                        n_distance_cost_ij = n_distance_cost_ij + len(distance_cost_ij)

                                                        sum_of_duration_cost_ij = sum_of_duration_cost_ij + np.sum(duration_cost_ij)
                                                        n_duration_cost_ij = n_duration_cost_ij + len(duration_cost_ij)

                                                        sum_of_ground_damage_cost_ij = sum_of_ground_damage_cost_ij + np.sum(ground_damage_cost_ij)
                                                        n_ground_damage_cost_ij = n_ground_damage_cost_ij + len(ground_damage_cost_ij)

                                                else:

                                                        #
                                                        # If candidate list is depleted, choose the next node to be the one
                                                        # with the largest desirability eta_ij
                                                        #

                                                        distance_cost_ij = distance_matrix[current_node, unvisited_pickup_nodes] / estimate_of_distance_cost_in_heuristic_desirability
                                                        duration_cost_ij = distance_matrix[current_node, unvisited_pickup_nodes] / parameters.forwarder_speed / estimate_of_duration_cost_in_heuristic_desirability
                                                        ground_damage_cost_ij = ground_damage_matrix[current_node, unvisited_pickup_nodes]*(1.0 + forwarder_load_volume / parameters.ground_damage_model_v0) / estimate_of_ground_damage_cost_in_heuristic_desirability

                                                        total_cost_ij = parameters.w_length*distance_cost_ij + parameters.w_duration*duration_cost_ij + parameters.w_ground_damage*ground_damage_cost_ij                                                        
                                                        eta_ij = 1.0 / total_cost_ij

                                                        i_next_node = np.argmax(eta_ij)
                                                        next_node = unvisited_pickup_nodes[i_next_node]


                                                        #
                                                        # Update sums for the components of desirability
                                                        #

                                                        sum_of_distance_cost_ij = sum_of_distance_cost_ij + np.sum(distance_cost_ij)
                                                        n_distance_cost_ij = n_distance_cost_ij + len(distance_cost_ij)

                                                        sum_of_duration_cost_ij = sum_of_duration_cost_ij + np.sum(duration_cost_ij)
                                                        n_duration_cost_ij = n_duration_cost_ij + len(duration_cost_ij)

                                                        sum_of_ground_damage_cost_ij = sum_of_ground_damage_cost_ij + np.sum(ground_damage_cost_ij)
                                                        n_ground_damage_cost_ij = n_ground_damage_cost_ij + len(ground_damage_cost_ij)


                                        #
                                        # Add the new node to the path and remove it from the list of unvisited nodes
                                        #

                                        this_path.append(next_node)
                                        unvisited_pickup_nodes.remove(next_node)

                                        #
                                        # Update pheromone concentration for this edge
                                        #

                                        pheromone_matrix[current_node, next_node] = (1.0 - parameters.rho) * pheromone_matrix[current_node, next_node] + parameters.rho * this_subproblem_tau_0


                                        #
                                        # Update forwarder load volume
                                        #

                                        forwarder_load_volume = forwarder_load_volume + wood_matrix[next_node, 1]

                                        current_node = next_node

                                #
                                # The forwarder has run out of capacity, or all nodes have been visited, so we must return to the roadside node to unload
                                #

                                next_node = 0
                                this_path.append(next_node)

                                #
                                # Update pheromone concentration for this edge
                                #

                                pheromone_matrix[current_node, 0] = (1.0 - parameters.rho) * pheromone_matrix[current_node, next_node] + parameters.rho * this_subproblem_tau_0

                                current_node = next_node


                        #
                        # A new solution has been created. Compute its cost, length, and duration.
                        #

                        t_begin_solution = time.time()

                        this_solution = smart_log_classes_aco.solution(this_path, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf, np.Inf)
                        this_solution.weighted_cost, this_solution.length, this_solution.duration, this_solution.ground_damage, this_solution.standardized_length, this_solution.standardized_duration, this_solution.standardized_ground_damage, this_solution.n_types_of_wood_by_load = functions.compute_cost_of_solution(this_solution, distance_matrix, ground_damage_matrix, wood_matrix, estimated_mean_length_of_path, estimated_mean_duration_of_path, estimated_mean_ground_damage_of_path)

                        t_end_solution = time.time()
                        time_for_computing_solution_cost_for_this_subproblem = time_for_computing_solution_cost_for_this_subproblem + (t_end_solution - t_begin_solution)

                        print("Iteration number %d, ant number = %d, weighted total cost = w_length*%f + w_duration*%f + w_ground_damage*%f = %f, length = %f m, duration = %f s, damage = %f" % (i_iteration + 1, i_ant + 1, this_solution.standardized_length, this_solution.standardized_duration, this_solution.standardized_ground_damage, this_solution.weighted_cost, this_solution.length, this_solution.duration, this_solution.ground_damage))


                        #
                        # See if we found a new champion path
                        #

                        if this_solution.weighted_cost < best_solution_so_far.weighted_cost:

                                if best_solution_so_far.path == []:

                                        print("---> New champion found at iteration %d, ant %d! Cost is %f" % (i_iteration + 1, i_ant + 1, this_solution.weighted_cost))

                                else:

                                        print("---> New champion found at iteration %d, ant %d! Old minimum cost was %f, new one is %f" % (i_iteration + 1, i_ant + 1, best_solution_so_far.weighted_cost, this_solution.weighted_cost))

                                best_solution_so_far = this_solution


                #
                # After all ants have completed a tour, update the pheromone
                # concentrations using the best solution so far
                #

                for i in range(0, len(best_solution_so_far.path) - 1):

                        this_node = best_solution_so_far.path[i]
                        next_node = best_solution_so_far.path[i+1]
                        pheromone_matrix[this_node, next_node] = (1.0 - parameters.rho) * pheromone_matrix[this_node, next_node] + parameters.rho / best_solution_so_far.weighted_cost


                #
                # Keep track of the best solution cost vs. iteration
                #

                best_solution_cost_vs_iteration[i_iteration, 0] = i_iteration + 1
                best_solution_cost_vs_iteration[i_iteration, 1] = best_solution_so_far.weighted_cost


        #
        # The optimization of this subproblem is now completed
        #

        t_end_optimization_of_this_subproblem = time.time()
        t_end_subproblem = time.time()

        time_for_optimization_of_this_subproblem = t_end_optimization_of_this_subproblem - t_begin_optimization_of_this_subproblem
        total_time_for_optimization = total_time_for_optimization + time_for_optimization_of_this_subproblem

        time_for_this_subproblem = t_end_subproblem - t_begin_subproblem

        #
        # Compute the mean number of piles in a tour of the best solution of
        # this subproblem, and update the total counters for this
        #

        best_solution_path = best_solution_so_far.path
        number_of_tours_in_best_solution = best_solution_path.count(0) - 1
        number_of_piles_in_best_solution = len(best_solution_path) - best_solution_path.count(0)
        mean_number_of_piles_per_tour_in_best_solution = number_of_piles_in_best_solution / number_of_tours_in_best_solution                

        total_number_of_tours = total_number_of_tours + number_of_tours_in_best_solution
        total_number_of_piles = total_number_of_piles + number_of_piles_in_best_solution

        #
        # Compute also the total m**3 forwarded in the best solution
        # of this subproblem, and update the total counter for this
        #

        volume_forwarded_in_this_subproblem = 0.0
        
        for i_node in best_solution_so_far.path:
        
                volume_forwarded_in_this_subproblem = volume_forwarded_in_this_subproblem + wood_matrix[i_node, 1]
        
        total_volume_of_wood_forwarded = total_volume_of_wood_forwarded + volume_forwarded_in_this_subproblem
        
        print("")
        print("===> Optimization completed for subproblem %d <===" % i_subproblem)
        print("")
        print("The best solution is the following path: ", best_solution_so_far.path)
        print("")
        print("The number of unique wood types for each load in this solution is as follows: ", best_solution_so_far.n_types_of_wood_by_load)
        print("")
        print("This solution has a total weighted cost of w_length*%f + w_duration*%f + w_ground_damage*%f = %f, a length of %f m or %f km, a duration of %f s or %f h, and total ground damage of %f" % (best_solution_so_far.standardized_length, best_solution_so_far.standardized_duration, best_solution_so_far.standardized_ground_damage, best_solution_so_far.weighted_cost, best_solution_so_far.length, best_solution_so_far.length / 1000.0, best_solution_so_far.duration, best_solution_so_far.duration / 3600.0, best_solution_so_far.ground_damage))
        print("")
        print("Mean number of piles per tour in this best solution was %f" % mean_number_of_piles_per_tour_in_best_solution)
        print("")
        print("Volume of wood forwarded in this subproblem was %f m**3" % volume_forwarded_in_this_subproblem)
        print("")

        print("Wall-clock time consumed by this subproblem:")
        print("")
        print("Total %f s, of which:" % time_for_this_subproblem)
        print("")
        print("---> Computing static cost matrices %f s" % time_for_computing_cost_matrices_for_this_subproblem)
        print("---> Optimization loop %f s" % time_for_optimization_of_this_subproblem)
        print("     ---> of which, computing solution costs %f s" % time_for_computing_solution_cost_for_this_subproblem)

        #
        # Update the total cost, length, and duration of the solution
        # to the full problem
        #

        total_weighted_cost = total_weighted_cost + best_solution_so_far.weighted_cost
        total_length = total_length + best_solution_so_far.length
        total_duration = total_duration + best_solution_so_far.duration
        total_ground_damage = total_ground_damage + best_solution_so_far.ground_damage

        
        #
        # Plot results for this subproblem
        #
        
        #
        # Plot the best solution on the real-life graph
        #

        functions.plot_best_solution(the_real_life_graph, best_solution_so_far, i_subproblem, pile_node_colors, pile_node_sizes, pile_positions_dictionary, pile_color_legend_handles, graph_axis_limits, subproblem)

        #
        # Plot best solution cost vs. iteration
        #

        functions.plot_best_solution_cost_vs_iteration(i_subproblem, best_solution_cost_vs_iteration)

        #
        # Plot the load composition vs. tour number in the best solution
        #

        functions.plot_load_composition_vs_tour(best_solution_so_far, wood_matrix, wood_type_dictionary, new_wood_type_dictionary, new_wood_type_color_dictionary, i_subproblem)


        #
        # Compute and print out statistics on the components of the heuristic
        # desirability
        #

        print("")
        print("Statistics on components of eta_ij:")
        print("")
        print("D_ij: mean = %f" % (sum_of_distance_cost_ij / n_distance_cost_ij))
        print("T_ij: mean = %f" % (sum_of_duration_cost_ij / n_duration_cost_ij))
        print("S_ij: mean = %f" % (sum_of_ground_damage_cost_ij / n_ground_damage_cost_ij))
        
        
        i_subproblem = i_subproblem + 1

        
#
# Loop over subproblems is complete
#

print("")
print("*** Loop over subproblems is complete ***")

print("")
print("The total solution for the full problem has a total weighted cost of %f, a length of %f m or %f km, a duration of %f s or %f h, and total ground damage of %f" % (total_weighted_cost, total_length, total_length / 1000.0, total_duration, total_duration / 3600.0, total_ground_damage))
print("")
print("The total volume of forwarded wood was %f m**3" % total_volume_of_wood_forwarded)
print("")
print("The total number of piles was %d, the total number of tours was %d, and the mean number of piles per tour was %f" % (total_number_of_piles, total_number_of_tours, total_number_of_piles / total_number_of_tours))


#
# Compute and report final timings for the entire program
#

t_end_program = time.time()
time_total_for_program = t_end_program - t_begin_program

print("")
print("Wall-clock time consumed by the entire program:")
print("")
print("Total %f s, of which:" % time_total_for_program)
print("")
print("---> Creating real-life graph %f s" % time_for_creating_real_life_graph)
print("---> Creating decomposition of real-life graph %f s" % time_for_creating_decomposition)
print("---> Computing cost matrices %f s" % total_time_for_computing_cost_matrices)
print("---> Optimizing loop %f s" % total_time_for_optimization)


#
# All done. Exit.
#

print("")
print("All done. Exiting.")
print("")


exit(0)
