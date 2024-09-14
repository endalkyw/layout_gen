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


# generate_cm_layouts(10, 10, 2, "N", source_first=True)

lib = gdstk.Library(name="LibraryName", unit=1e-9, precision=1e-12)
cell_i = lib.new_cell("base_cell")

# g_t = "cm"
g_t = "dp_0"
# g_t = "dp_1"
g_t = "dp_2"

create_base_layout(cell_i, 12, 6, "N", gate_type=g_t, stack=1, multiplier=4, orientation="V", body_contact = False)
write_gds(cell_i, "base_cell")


# For single mos test
# m0 = Mos({'id': 'A', 'fins': 10, 'fingers': 19, 'stack': 2, 'multiplier': 3, 'mos_type': "N"})
# cell = create_mos(m0, labels=["d", 0, 0, "s"], con=[1,1,1,1], orientation="V")
# write_gds(cell, "nmos")


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