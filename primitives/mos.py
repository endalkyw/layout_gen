from utilities.elements_utils import *
from utilities.pvector import P
from core import create_base_layout, create_fabric
from core.constants import Const as C
from utilities.via import *
import math

def create_mos(m,
               fabric_on = True,
               con= [2, 1, 2, 1],
               labels = [0, 0, 0, 0],
               orientation = "v",
               gate_connection_type = "cm"):

    cell_i = create_empty_cell("single_nmos", unit=1e-9,  precision=1e-12)
    fin_dummy  = lp['fins']['dummies']
    fing_dummy = lp['poly']['dummies']
    fin_no = m.fins + 4 * fin_dummy - 2
    poly_no = m.stack*m.fingers + 2 * fing_dummy
    poly_height = fin_no * lp['fins']['pitch']
    fin_len = (poly_no - 1) * lp['poly']['pitch']
    points = create_base_layout(cell_i, m.fingers, m.fins, m.mos_type, stack=m.stack, multiplier = m.multiplier, orientation=orientation)

    for pnt in points:
        loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
        n0 = m.fingers//2 + int(m.fingers % 2)

        # metal connections
        src_drn_cons = con[0] + con[2]
        avail_v_space = loc["y_span"][1] - loc["y_span"][0]
        sd_space = (src_drn_cons-1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2*C.m_m_ext
        v_sp  = (avail_v_space - sd_space)
        v_sp2 = avail_v_space - ((con[2]-1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2*C.m_m_ext)
        # drain metal and via array

        if v_sp < 0:
            add_v_metal_array(cell_i, P(C.rx_m1 + m.stack * lp['M1']['pitch'], v_sp), abs(v_sp - 1), m.stack * 2 * lp['M2']['pitch'], n0, "M1")
            if v_sp2 < 0:
                add_v_metal_array(cell_i, P(C.rx_m1, v_sp2), abs(v_sp2 - 1),
                                  m.stack * 2 * lp['M2']['pitch'], n0+(m.fingers+1) % 2, "M1")
        else:
            v_sp /= 2
        for i in range(con[0]):
            add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1 + m.stack*lp['M1']['pitch'], pnt["p0"][0][1] + v_sp + i*lp["M2"]['pitch']+ C.m_m_ext),  m.stack*2*lp['M2']['pitch'], n0, 'H', "V1")
            add_metal(cell_i, P(pnt["p0"][0][0] + C.rx_m1 + m.stack*lp['M1']['pitch'] - C.m_v_ext, pnt["p0"][0][1] + v_sp + i*lp["M2"]['pitch']+ C.m_m_ext), m.stack*(n0-1)*2*lp['M2']['pitch']+lp['M1']['width']+2*C.m_v_ext, "H", "M2")

        # source metal and via array
        for i in range(con[2]):
            add_metal(cell_i, P(loc["x_span"][0],  pnt["p0"][0][1] + v_sp + con[0]*lp["M2"]['pitch'] + i*lp["M2"]['pitch']+ C.m_m_ext), loc["x_span"][1]-loc["x_span"][0], "H", "M2")
            add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1, pnt["p0"][0][1] + v_sp + con[0]*lp["M2"]['pitch'] + i*lp["M2"]['pitch']+ C.m_m_ext), m.stack*2*lp['M1']['pitch'], m.fingers//2+1, 'H', "V1")

        # gate metal and via array
        add_metal(cell_i, P(loc["x_span"][0]+C.rx_m1, pnt["g"][0][1]), loc["x_span"][1]-loc["x_span"][0]-2*C.rx_m1, "H", "M1")
        # g.update(0, loc["g"].y)




    p = cell_i.bounding_box()
    poly_n = math.ceil((p[1][0] - p[0][0])/lp["poly"]["pitch"])-1
    fin_n = math.ceil((p[1][1] - p[0][1])//lp["fins"]["pitch"])
    print(poly_n, fin_n)
    x_ = 0
    y_ = 0
    cell = create_empty_cell("cm", 1e-9, 1e-12)
    if fabric_on:
        create_fabric(cell, poly_n, fin_n)
        x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
        y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2
    add_transformed_polygons(cell_i, cell, (x_, y_))


    return cell






