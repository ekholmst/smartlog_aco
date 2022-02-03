# smartlog_aco
SmartLog-ACO is a heuristic, multi-objective approach to solving the forwarder routing problem in cut-to-length forest harvesting operations. See E. Holmström, J. Nikander, J. Backman, K. Väätäinen, J. Uusitalo, and P. Jylhä, "A multi-objective optimization strategy for timber forwarding in cut-to-length harvesting operations" (in preparation) for details.

Use the following command to run the code:

smart_log_aco.py pile_data.txt path_data.txt

where pile_data.txt is the path to a file containing the data on log piles in the following format:

<northing (m)> <easting (m)> <assortment (an integer value)> <total volume of the pile (m**3)> <number of logs in the pile> <maximum distance between two logs in the pile (m)>
  
and path_data.txt is the path to a file containing the forest machine GNSS trace, used in constructing the network of trails and piles, and has the following format:
  
<altitude (m)> <northing (m)> <easting (m)>

Uses the following packages: sys, os, string, math, re, numpy, matplotlib, math, warnings, itertools, copy, random, networkx, scipy, PIL, time.

Contact: eero.holmstrom@luke.fi
