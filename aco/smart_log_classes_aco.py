#
# -*- coding: utf-8 -*-
#
# smart_log_classes_aco.py
#
# User-defined classes for smart_log_aco.py.
#
# Copyright 2022 Natural Resources Institute Finland and Eero
# Holmström (eero.holmstrom@luke.fi)
#
# This program is distributed under the terms of the GNU Lesser
# General Public License version 3.0 or any later version.
#


#
# A class describing a pick-up node, i.e., a pile of wood to be picked up and
# transported to a roadside position.
#

class smart_log_pickup_node:
        
        def __init__(self, x_position, y_position, type_of_wood, amount_of_wood):
                
                self.x_position = x_position
                self.y_position = y_position
                self.wood_pile = smart_log_wood_pile(type_of_wood, amount_of_wood)
                


#
# A class describing a pile of wood (on the ground or already loaded onto the forwarder).
#

class smart_log_wood_pile:

        def __init__(self, type_of_wood, amount_of_wood):

                self.type_of_wood = type_of_wood
                self.amount_of_wood = amount_of_wood

                

#
# A class describing a roadside delivery node, i.e., where the piles of wood
# are to be brought to.
#

class smart_log_roadside_node:
        
        def __init__(self, x_position, y_position):
                
                self.x_position = x_position
                self.y_position = y_position


#
# A class describing a solution, which consists of a path through the nodes and
# the corresponding loadings or unloadings of wood.
#

class solution:
                
        def __init__(self, path, weighted_cost, length, duration, ground_damage, standardized_length, standardized_duration, standardized_ground_damage):

                # The sequence of nodes to go through
                self.path = path

                # Total weighted cost of this path
                self.weighted_cost = weighted_cost

                # Total length of this path (in m)
                self.length = length

                # Total time required to complete this path (in s)
                self.duration = duration

                # Total ground damage that completing this path causes
                self.ground_damage = ground_damage

                # Total standardized length of this path
                self.standardized_length = standardized_length

                # Total standardized time required to complete this path
                self.standardized_duration = standardized_duration

                # Total standardized ground damage that completing this path causes
                self.standardized_ground_damage = standardized_ground_damage
                
                # The number of unique types of wood for each load in this solution
                self.n_types_of_wood_by_load = []
                
