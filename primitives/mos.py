from utilities.elements_utils import *
from utilities.pvector import P
from core import create_base_layout, create_fabric
from core.constants import Const as C
from utilities.via import *

def create_mos(m,
               fabric_on = True,
               con= [2, 1, 2, 1],
               labels = [0, 0, 0, 0],
               orientation = "v",
               gate_connection_type = "cm"):

    def create_cell(bulk_on):
        cell_i = create_empty_cell("single_nmos", unit=1e-9,  precision=1e-12)
        fin_dummy  = lp['fins']['dummies']
        fing_dummy = lp['poly']['dummies']
        fin_no = m.fins + 4 * fin_dummy - 2
        poly_no = m.stack*m.fingers + 2 * fing_dummy
        poly_height = fin_no * lp['fins']['pitch']
        fin_len = (poly_no - 1) * lp['poly']['pitch']
        points = create_base_layout(cell_i, m.fingers, m.fins, m.mos_type, stack=m.stack, body_contact=bulk_on, gate_type=gate_connection_type)


        print(points["p0"])
        # loc = {"x_span": [points["p0"][0][0], points["p1"][0][0]], "y_span": [points["p0"][0][1], points["p1"][0][1]]}
        # # show_layout(cell_i,(2,1))
        # n0 = m.fingers//2 + int(m.fingers % 2)
        # # d = P(0, 0)
        # # s = P(0, 0)
        # # g = P(0, 0)
        # # b = P(0, 0)
        #
        # # b.update(loc["b"].x, loc["b"].y)
        # # metal connections
        # src_drn_cons = con[0] + con[2]
        # avail_v_space = loc["y_span"][1] - loc["y_span"][0]
        # sd_space = (src_drn_cons-1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2*C.m_m_ext
        # v_sp  = (avail_v_space - sd_space)
        # v_sp2 = avail_v_space - ((con[2]-1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2*C.m_m_ext)
        # # drain metal and via array
        # if v_sp < 0:
        #     add_v_metal_array(cell_i, P(C.rx_m1 + m.stack * lp['M1']['pitch'], v_sp), abs(v_sp - 1), m.stack * 2 * lp['M2']['pitch'], n0, "M1")
        #     if v_sp2 < 0:
        #         add_v_metal_array(cell_i, P(C.rx_m1, v_sp2), abs(v_sp2 - 1),
        #                           m.stack * 2 * lp['M2']['pitch'], n0+(m.fingers+1)%2, "M1")
        # else:
        #     v_sp /= 2
        # for i in range(con[0]):
        #     add_via_array(cell_i, (C.rx_m1 + m.stack*lp['M1']['pitch'], v_sp + i*lp["M2"]['pitch']+ C.m_m_ext),  m.stack*2*lp['M2']['pitch'], n0, 'H', "V1")
        #     add_metal(cell_i, P(C.rx_m1 + m.stack*lp['M1']['pitch'] - C.m_v_ext, v_sp + i*lp["M2"]['pitch']+ C.m_m_ext), m.stack*(n0-1)*2*lp['M2']['pitch']+lp['M1']['width']+2*C.m_v_ext, "H", "M2")
        #     # d.update(C.rx_m1 + m.stack*lp['M1']['pitch'] - C.m_v_ext, v_sp + i*lp["M2"]['pitch']+ C.m_m_ext)
        # # source metal and via array
        # for i in range(con[2]):
        #     add_metal(cell_i, P(loc["x_span"][0],  v_sp + con[0]*lp["M2"]['pitch'] + i*lp["M2"]['pitch']+ C.m_m_ext), loc["x_span"][1]-loc["x_span"][0], "H", "M2")
        #     add_via_array(cell_i, (C.rx_m1, v_sp + con[0]*lp["M2"]['pitch'] + i*lp["M2"]['pitch']+ C.m_m_ext), m.stack*2*lp['M1']['pitch'], m.fingers//2+1, 'H', "V1")
        #     # s.update(loc["x_span"][0],  v_sp + con[0]*lp["M2"]['pitch'] + i*lp["M2"]['pitch']+ C.m_m_ext)
        #
        # # gate metal and via array
        # add_metal(cell_i, P(loc["x_span"][0]+C.rx_m1, loc["g"].y),loc["x_span"][1]-loc["x_span"][0]-2*C.rx_m1, "H", "M1")
        # # g.update(0, loc["g"].y)
        return cell_i




    cell = create_empty_cell("single_nmos", unit=1e-9, precision=1e-12)
    delta = 0
    for i in range(m.multiplier):
        cell_i = create_cell(bulk_on = (i + 1)//m.multiplier)
        p = cell_i.bounding_box()
        print(p)
        if orientation == "V":
            add_transformed_polygons(cell_i, cell, (0, delta))
            delta += p[1][1]

        elif orientation == "H":
            add_transformed_polygons(cell_i, cell, (delta, 0))
            delta += p[1][0]



    # x_ = 0
    # y_ = 0
    #
    #
    # cell = create_empty_cell("cm", 1e-9, 1e-12)
    # if fabric_on:
    #     create_fabric(cell, poly_no, fin_no)
    #     x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
    #     y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2
    #
    # add_transformed_polygons(cell_i, cell, (x_, y_))
    # # show_layout(cell,(2,2))
    # # write_gds(cell, "end"+str(1))
    #
    #
    # d.transform_p(x_, y_)
    # s.transform_p(x_, y_)
    # g.transform_p(x_+16, y_+16)
    # b.transform_p(x_+16, y_+16)
    #
    # if labels[0] != 0:
    #     add_label(cell, labels[0], d, "M1LABEL")
    #
    # if labels[1] != 0:
    #     add_label(cell, labels[1], g, "M1LABEL")
    #
    # if labels[2] != 0:
    #     add_label(cell, labels[2], s, "M1LABEL")
    #
    # if labels[3] != 0:
    #     add_label(cell, labels[3], b, "M1LABEL")

    return cell






