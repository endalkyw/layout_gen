import numpy as np
from joblib import Memory
from primitives.mos import *
import gdstk
from primitives.differential_pair import differential_pair
from utilities.elements_utils import *
from primitives.current_mirror import current_mirror
from core import create_base_layout, create_fabric
import math
import time
from functools import wraps
from multiple_layout_generator import *
from core.base_cell import *


# # ----------------------- create base layout test---------------------------------------------
# cell = create_empty_cell("base_cell",unit=1e-9, precision=1e-12)
# inputs = {
#     "cell"          :   cell,
#     "total_fingers" :   12,
#     "total_fins"    :   10,
#     "mos_type"      :   "N",
#     "stack"         :   1,
#     "multiplier"    :   10,
#     "orientation"   :   "VH",
#     "gate_type"     :   "cm",    # "cm" "dp_0" "dp_1" "dp_2"
#     "body_contact"  :   False
# }
# points = create_base_layout(**inputs)
# write_gds(cell, "base_cell")
#
# for i in points:
#     print(i)
## ---------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# For single mos test
m0 = Mos({'id': 'A', 'fins': 10, 'fingers': 18, 'stack': 1, 'multiplier': 5, 'mos_type': "N"})
cell = create_mos(m0, labels=["d", 0, 0, "s"], con=[2, 2, 2, 1], orientation="V", fabric_on=True)
write_gds(cell, "nmos")


# plt.rcParams['figure.facecolor'] = 'black'  # Set figure background
# plt.rcParams['axes.facecolor'] = 'black'    # Set axes background
# show_layout(cell, fig_size=(6,4))
# ---------------------------------------------------------------------------------------------





# For single current mirror test
# m0 = Mos({'id': 'A', 'fins': 10, 'fingers': 4, 'stack': 1, 'multiplier': 1, 'mos_type': "N"})
# m1 = Mos({'id': 'A', 'fins': 10, 'fingers': 4, 'stack': 1, 'multiplier': 1, 'mos_type': "N"})
# c1 = current_mirror(m0, m1, "test_cm")
# c1.create_layout(1, ("d0", "d1", "s"), [1, 1, 1], 1, True)
# write_gds(c1.cell, "nmos")


# For single differential pair test
# m0 = Mos({'id': 'A', 'fins': 10, 'fingers': 4, 'stack': 2, 'multiplier': 1, 'mos_type': "N"})
# m1 = Mos({'id': 'A', 'fins': 10, 'fingers': 4, 'stack': 2, 'multiplier': 1, 'mos_type': "N"})
# c1 = differential_pair(m0, m1, "test_dp")
# c1.create_layout(0)
# write_gds(c1.cell, "nmos_dp")




# generate_cm_layouts(10, 10, 2, "N", source_first=True)
