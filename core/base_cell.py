from utilities.elements_utils import *
from utilities.via import *
from utilities.pvector import P
import gdstk
from core.constants import Const as C
from core.base_fabric import ConnectionRects as CR


def transform_contact_points(lists, x_offset, y_offset):
    for key, item in lists.items():
       if key in ["b", "s", "d", "g", "p0", "p1"]:
        for ls in item:
          ls[0] += x_offset
          ls[1] += y_offset

    # print(lists)

def create_base_layout(cell, total_fingers, total_fins, mos_type, stack = 1, multiplier = 1, orientation = "V", gate_type="cm", body_contact = True):
    x_i = 0
    y_i = 0
    x_f = (stack*total_fingers + 1) * lp['poly']['pitch']
    y_f = total_fins * lp['fins']['pitch']

    def create_cell(cell_i, bulk_on):
        # add active region
        contact_points = {"b": [], "s": [], "d": [], "g": [], "p0": [], "p1": [], "m": multiplier, "o": orientation}

        contact_points["p0"].append([x_i, y_i])
        contact_points["p1"].append([x_f, y_f])

        add_region(cell_i, (x_i, y_i), (x_f, y_f), "RX")
        add_region(cell_i, (x_i, y_i), (x_f, y_f), "UNK_1")  # 212:121
        add_region(cell_i, (x_i, y_i), (x_f, y_f), "TT")

        h_ext = (lp["poly"]["pitch"] - lp["poly"]["width"]) / 2 + lp["poly"]["width"] / 2
        v_ext = lp["fins"]["pitch"]

        add_region(cell_i, (x_i - h_ext, y_i - v_ext), (x_f + h_ext, y_f + v_ext), "RVT")

        if mos_type == "P":
            add_region(cell_i, (x_i - h_ext - 9, y_i - v_ext - 4), (x_f + h_ext + 9, y_f + v_ext + 4), "BP")

        add_region(cell_i, (x_i, y_i), (x_f, y_f), "UNK")
        add_region(cell_i, (x_i - lp['poly']['width'] / 2, y_i - lp["AUXPC"]["ext_len"]),
                   (x_i + lp['poly']['width'] / 2, y_f + lp["AUXPC"]["ext_len"]), "AUXPC")
        add_region(cell_i, (x_f - lp['poly']['width'] / 2, y_i - lp["AUXPC"]["ext_len"]),
                   (x_f + lp['poly']['width'] / 2, y_f + lp["AUXPC"]["ext_len"]), "AUXPC")


        # Adding body contact
        if bulk_on:
            b_x = (x_i + x_f) / 2
            if gate_type == "cm" or gate_type == "dp_0":
                b_y = y_f + C.bulk_h1
            else:
                b_y = y_f + C.bulk_h2

            bulk_w = max(150, 0.9 * (stack*total_fingers * lp['poly']['pitch']))
            bulk_h = 100

            m1 = 20
            m3 = 25
            m2 = 30
            m4 = 35

            if mos_type == "N":
                add_region(cell_i, (b_x - bulk_w / 2, b_y), (b_x + bulk_w / 2, b_y + bulk_h), "BP")

            add_region(cell_i, (b_x - bulk_w / 2 + m1, b_y + m1), (b_x + bulk_w / 2 - m1, b_y + bulk_h - m1), "RX")
            add_region(cell_i, (b_x - bulk_w / 2 + m2, b_y + m3), (b_x + bulk_w / 2 - m2, b_y + bulk_h - m3), "TT")
            add_region(cell_i, (b_x - bulk_w / 2 + m4, b_y + m2), (b_x + bulk_w / 2 - m4, b_y + bulk_h - m2), "CA")

            r0 = [(b_x - bulk_w / 2 + m4, b_y + m2), (b_x + bulk_w / 2 - m4, b_y + bulk_h - m2)]
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["b"].append([r0[0][0], r0[0][1]])

        # Adding CA (14:0) and M1 (15:0) [Drain and Source diffusions]
        add_v_metal_array(cell_i, P(C.rx_m1-1,0), y_f, lp["CA"]["pitch"], stack*total_fingers+1, "CA")
        add_v_metal_array(cell_i, P(C.rx_m1,0), y_f,  stack*lp["M1"]["pitch"], total_fingers+1, "M1")
        fill_with_vias_array(cell_i, "V",(C.rx_m1,0),(C.rx_m1+lp["M1"]["width"],y_f),
                             total_fingers+1, stack*lp["CA"]["pitch"], "V0")

        contact_points["s"].append([C.rx_m1, y_i])
        contact_points["d"].append([C.rx_m1+lp["M1"]["pitch"], y_i])


        # add gate connectors
        # add_poly_m1_connector(cell_i, total_fingers, y_f + 65)
        x_i_ = lp["poly"]["pitch"] - lp["M1"]["width"] / 2
        pitch = lp["poly"]["pitch"]

        # gate cb (245/0) <- from the active layer (2/0)
        d1 = x_i + C.rx_m1 + lp["M1"]["width"]
        l_x = stack*lp["poly"]["pitch"]*total_fingers
        d2 = d1 + l_x / 2 - lp["M1"]["width"]
        d3 = d2 + lp["M1"]["width"]
        d4 = d3 + l_x / 2 - lp["M1"]["width"]

        if gate_type == "cm":
            r0 = (d1, y_f + C.gate_h0), (d1 + l_x - lp["M1"]["width"], y_f + C.gate_h0+lp["CB"]["width"])
            add_region(cell_i, r0[0], r0[1], "CB")
            fill_area_vias(cell_i, r0, "V0")
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["g"].append(list(r0[0]))

        elif gate_type == "dp_0": # clustered type 2 gates
            r0 = (d1, y_f + C.gate_h0), (d2 , y_f + C.gate_h0+lp["CB"]["width"])
            add_region(cell_i, r0[0], r0[1], "CB")
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["g"].append(list(r0[0]))

            r1 = (d3, y_f + C.gate_h0), (d4, y_f + C.gate_h0+lp["CB"]["width"]),
            add_region(cell_i,  r1[0], r1[1], "CB")
            fill_area_vias(cell_i, r1, "V0")
            fill_v_metals(cell_i, r1, "M1")
            contact_points["g"].append(list(r1[0]))

        elif gate_type == "dp_1":  # for interdigitated
            pitch = 2*pitch*stack

            r0 = [(d1, y_f + C.gate_h0), (d1 + stack*2*lp["M1"]["pitch"] - lp["M1"]["width"], y_f + C.gate_h0 + 30)]
            add_region_array(cell_i, r0[0], (r0[1][0]-r0[0][0], r0[1][1]-r0[0][1]), 2*pitch, max(1,total_fingers//4), "H", "CB")
            fill_with_vias_array(cell_i, "H", r0[0], r0[1], max(1,total_fingers//4), 2*pitch, "V0")
            fill_with_v_M1_array(cell_i, r0, 2*pitch,max(1,total_fingers//4), "M1")
            contact_points["g"].append(list(r0[0]))

            r1 = [(d1 + stack*2*lp["M1"]["pitch"], y_f + C.gate_h1), (d1 + stack*4*lp["M1"]["pitch"] - lp["M1"]["width"], y_f + C.gate_h1 + 30)]
            add_region_array(cell_i, r1[0], (r1[1][0]-r1[0][0], r1[1][1]-r1[0][1]), 2*pitch, max(1,total_fingers//4), "H", "CB")
            fill_with_vias_array(cell_i, "H", r1[0], r1[1], max(1,total_fingers//4), 2*pitch, "V0")
            fill_with_v_M1_array(cell_i, r1, 2*pitch,max(1,total_fingers//4), "M1")
            contact_points["g"].append(list(r1[0]))

        elif gate_type == "dp_2":  # for common centroid
            n1_ = total_fingers//4
            n2_ = total_fingers//2
            pitch *= stack

            r0 = [(d1, y_f + C.gate_h0),  (d1 + l_x // 4 - lp["M1"]["width"], y_f + C.gate_h0 + 30)]
            add_region(cell_i, r0[0], r0[1], "CB")
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["g"].append(list(r0[0]))

            k1 = stack*n1_ * lp['poly']['pitch']
            r1 = [(d1 + k1, y_f + C.gate_h1),  (d1 + k1 + l_x / 2 - lp["M1"]["width"], y_f + 30+ C.gate_h1)]
            add_region(cell_i, r1[0], r1[1], "CB")
            fill_area_vias(cell_i, r1, "V0")
            fill_v_metals(cell_i, r1, "M1")
            contact_points["g"].append(list(r1[0]))

            k2 = stack*(n1_ + n2_)* lp['poly']['pitch']
            r2 = [(d1 + k2, y_f + C.gate_h0),  (d1 + k2 + l_x / 4 - lp["M1"]["width"], y_f + C.gate_h0 + 30)]
            add_region(cell_i, r2[0], r2[1], "CB")
            fill_area_vias(cell_i, r2, "V0")
            fill_v_metals(cell_i, r2, "M1")
            contact_points["g"].append(list(r2[0]))

        elif gate_type == "dp_1x0":
            pitch *= stack
            r0 = [(x_i_-17.5, y_f + C.gate_h0), (x_i_ + pitch - 31, y_f +  C.gate_h0 + 30)]
            add_region(cell_i, r0[0], r0[1], "CB")
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["g"].append(list(r0[0]))

        elif gate_type == "dp_1x1":
            pitch *= stack
            r0 = [(x_i_-17.5, y_f + C.gate_h1), (x_i_ + pitch - 31, y_f +  C.gate_h1 + 30)]
            add_region(cell_i, r0[0], r0[1], "CB")
            fill_area_vias(cell_i, r0, "V0")
            fill_v_metals(cell_i, r0, "M1")
            contact_points["g"].append(list(r0[0]))

        return contact_points



    delta  = 0
    anchor_points =  []
    for i in range(multiplier):
        cell_i = create_empty_cell("single_nmos", unit=1e-9, precision=1e-12)
        ext_len = (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 46

        if orientation == "V":
            contact_points = create_cell(cell_i, bulk_on=(i + 1) // multiplier)
            p = cell_i.bounding_box()

            # bottom CT
            add_region(cell_i, (x_i - ext_len, y_i - C.ct_hb - lp['CT']['height']) , (x_f + ext_len, y_i - C.ct_hb), "CT")
            if i == multiplier-1:
                # top CT
                if gate_type == "cm" or gate_type == "dp_0":
                    add_region(cell_i, (x_i - ext_len, y_f + C.ct_h0), (x_f + ext_len, y_f + C.ct_h0 + lp['CT']['height']), "CT")
                else:
                    add_region(cell_i, (x_i - ext_len, y_f + C.ct_h1), (x_f + ext_len, y_f + C.ct_h1 + lp['CT']['height']), "CT")

            add_transformed_polygons(cell_i, cell, (0, delta))
            transform_contact_points(contact_points, 0, delta)
            anchor_points.append(contact_points)

            delta += p[1][1] + C.v_mult_space*lp["fins"]["pitch"]


        elif orientation == "H":
            contact_points = create_cell(cell_i, bulk_on=True)
            p = cell_i.bounding_box()

            # bottom CT
            add_region(cell_i, (x_i - ext_len, y_i - C.ct_hb - lp['CT']['height']), (x_f + ext_len, y_i - C.ct_hb),
                       "CT")
            # top CT
            if gate_type == "cm" or gate_type == "dp_0":
                add_region(cell_i, (x_i - ext_len, y_f + C.ct_h0),
                           (x_f + ext_len, y_f + C.ct_h0 + lp['CT']['height']), "CT")
            else:
                add_region(cell_i, (x_i - ext_len, y_f + C.ct_h1),
                               (x_f + ext_len, y_f + C.ct_h1 + lp['CT']['height']), "CT")\


            add_transformed_polygons(cell_i, cell, (delta, 0))
            transform_contact_points(contact_points, delta, 0)
            anchor_points.append(contact_points)

            delta += p[1][0] + C.h_mult_space*lp["poly"]["pitch"]


    p = cell.bounding_box()
    if mos_type == "P":
        add_region(cell, p[0], p[1], "NW")

    return anchor_points





