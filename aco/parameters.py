#
# -*- coding: utf-8 -*-
#
# parameters.py
#
# Parameters for smart_log_aco.py.
#
# Eero Holmström, 2019-2022
#


#
# General stuff
#

random_number_generator_seed = -1


#
# Harvestability map
#

harvestability_map_file = "/home/eero/workdir/smartlog/data/trafficability/whole_of_finland/KKL_SMK_Suomi_2020_10_26.tif"
harvestability_map_ulc_pixel_center_x_coordinate = 164008.0002122271
harvestability_map_ulc_pixel_center_y_coordinate = 7691992.0000002878
harvestability_map_pixel_size = 16.0


#
# Position of the roadside node in easting (m) and northing (m)
#

x_roadside_position = 
y_roadside_position = 


#
# Visualize an individual path
#

path_to_visualize = []
path_to_visualize_axis_limits = []


#
# Parameters for the forwarder
#

# Full loading capacity in m**3
forwarder_initial_capacity = 10.0

# Forwarder bunk is considered full when it has less than this amount of empty volume remaining
load_full_tolerance = 0.1

# Speed in m/s
forwarder_speed = 1.0


#
# Parameters for visualizations
# 

save_figures_to_files = True
plot_each_tour_separately = True
save_joining_and_merging_figures_to_files = False
save_individual_subgraph_plots_to_files = False

pickup_node_marker_scaling_factor = 2.0
default_pile_marker_size = 50.0
legend_wood_pile_marker_size = 35.0

roadside_node_scale = 100.0
road_point_markersize = 10.0
road_point_marker_edge_width = 3.0
raw_forwarder_data_markersize = 10.0
raw_forwarder_position_point_edge_width = 1.0

arrow_size = 25.0
path_width = 4.0
road_strip_width = 6.0
node_size_in_road_strip_graph = 150.0

fragment_path_width = 6.0
fragment_graph_node_size = 150.0

subgraph_path_width = 6.0
subgraph_cm_marker_size = 15.0
subgraph_cm_marker_edge_width = 3.0

node_size_in_before_and_after_joining_plots = 150.0

decomposition_alpha = 1.0

legend_region_graph_color_marker_size = 35.0

figure_size_h = 15.0
figure_size_v = 12.0
margin_for_pile_and_road_plots = 10

n_bins_wood_type_histogram = 20
pile_volume_histogram_bin_step = 0.25
pile_number_histogram_bin_step = 1.0
pile_extent_histogram_bin_step = 1.0
composition_bar_chart_width = 0.35

small_font_size = 10
default_font_size = 15
graph_node_font_size = 10
draw_graph_nodes_with_labels = False
title_pad = 25.0


#
# Parameters for the decomposition
#

maximum_number_of_nodes_per_region_graph = 50
minimum_number_of_nodes_per_subgraph = 1
decomposition_node_size = 65


#
# Smoothening of the input road data
#

road_point_identification_window_size = 1.0
road_point_identification_min_n_points_per_window = 1
n_road_points_for_pile_placement = 3
road_spline_parameter_step = 0.001


#
# Pile connectivity parameters. Angles are in degrees, distance in m.
#

pile_connectivity_initial_sector_range = 50.0
pile_connectivity_sector_range_increase = 10.0
pile_connectivity_max_distance_between_neighbors = 40.0
pile_connectivity_max_sector = 180.0


#
# Parameters for the time consumption of loading and unloading
#

#
# These are from Manner et al. (2013), Eq. (3)
#
t_load_a = 0.008336
t_load_b = 0.06435
t_load_c = 0.8678

#
# These are from Manner et al. (2013), Eq. (4)
#
t_unload_a = 0.00792
t_unload_b = 0.2312
t_unload_c = 0.0909


#
# Other parameters related to cost
#

minimum_cost_in_distance_matrix = 0.5
ground_damage_model_v0 = 10.0
cap_on_effective_number_of_unique_wood_types = 6


#
# Component weights in the optimization and in the heuristic desirability
#

w_length = 0.1
w_duration = 0.2
w_ground_damage = 0.3


#
# Parameters of the ACO algorithm
#

tau_0 = -1
rho = 0.01
alpha = 4
beta = 2
number_of_ants = 5
number_of_candidates = 20
max_iterations = 10
