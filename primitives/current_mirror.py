import gdstk

from utilities.elements_utils import *
from utilities.pvector import P
from core import create_base_layout, create_fabric
from primitives.mos import *
from utilities.via import *
import math
from core.constants import Const as C
from shapely.geometry import Polygon
import numpy as np

class current_mirror:
    def __init__(self, m0: Mos, m1: dict, cell_name: str):
        assert m0.fins >= 2, "Number of fins must be greater than 1"
        assert m0.stack == m1.stack
        assert m0.fins == m1.fins, "fins must be the same"
        assert m0.mos_type == m1.mos_type, "MOS type must be the same"

        self.m0 = m0
        self.m1 = m1
        self.total_fins = m0.fins
        self.stack = m0.stack
        self.total_fingers = m0.fingers + m1.fingers
        self.d0 = P(0, 0)
        self.d1 = P(0, 0)
        self.s  = P(0, 0)
        self.mos_type = m0.mos_type
        self.pin_positions = {"d0": [], "d1": [], "s": []}

        self.lib = gdstk.Library(name=cell_name, unit=1e-9, precision=1e-12)
        self.cell = self.lib.new_cell(cell_name)

    def create_layout(self, pattern: int,
                      labels = ["do", "d1", "s"], label_pos = ["c","c","c"],
                      con = [1, 1, 1], fabric_on=True):

        assert pattern == 0 or 1 or 2 or 3 or 4 or 5
        cell_i = create_empty_cell("cm", 1e-9, 1e-12)

        nd0 = self.m0.fingers // 2 + self.m0.fingers % 2
        nd1 = self.m1.fingers // 2 + self.m0.fingers % 2

        if pattern == 0 or pattern == 2 or pattern == 4:        # no diffusion break
            # common to all connections: source, gate, and the bulk is common (no need to add connection)
            points = create_base_layout(cell_i, self.total_fingers, self.total_fins, self.m0.mos_type,
                  stack=self.m0.stack, multiplier=self.m0.multiplier, orientation="V", body_contact=False)
            contact_rects = []

            for pnt in points:
                loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
                # metal connections
                if pattern == 0:
                    src_drn_cons = con[0] + con[2]
                    drn = con[0]
                else:
                    src_drn_cons = con[0] + con[1] + con[2]
                    drn = con[0]+con[1]

                avail_v_space = loc["y_span"][1] - loc["y_span"][0]
                sd_space = (src_drn_cons - 1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2 * C.m_m_ext
                v_sp = (avail_v_space - sd_space)/2

                # source metal and via array
                for i in range(con[2]):
                    xi, yi = loc["x_span"][0], pnt["p0"][0][1] + v_sp + drn * lp["M2"]['pitch'] + i * lp["M2"]['pitch'] + C.m_m_ext
                    rect_s = [[xi, yi], [xi + loc["x_span"][1] - loc["x_span"][0], yi + lp["M2"]['width']]]
                    add_metal(cell_i, P(xi, yi), loc["x_span"][1] - loc["x_span"][0], "H", "M2")
                    add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1,
                                           pnt["p0"][0][1] + v_sp + drn * lp["M2"]['pitch'] + i * lp["M2"][
                                               'pitch'] + C.m_m_ext), self.m0.stack * 2 * lp['M1']['pitch'],
                                  self.total_fingers // 2 + 1, 'H', "V1")
                    contact_rects.append(CR("s", rect_s))

                # gate connection
                rect_g = [(loc["x_span"][0] + C.rx_m1, pnt["g"][0][1] - 1), (
                loc["x_span"][0] + C.rx_m1 + loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1,
                pnt["g"][0][1] + lp["M1"]["width"] - 1)]
                add_vias_at_intersection(cell_i, rect_g, "M1")
                add_metal(cell_i, P(rect_g[0][0], rect_g[0][1]), loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1, "H","M2")
                contact_rects.append(CR("g", rect_g))


                # bulk metal and via array
                if pnt["b"]:
                    rect_b = [(loc["x_span"][0] + C.rx_m1, pnt["b"][0][1] + 4), (
                    loc["x_span"][0] + C.rx_m1 + loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1,
                    pnt["b"][0][1] + lp["M1"]["width"] + 4)]
                    add_vias_at_intersection(cell_i, rect_b, "M1")
                    add_metal(cell_i, P(rect_b[0][0], rect_b[0][1]), loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1, "H",
                              "M2")
                    contact_rects.append(CR("b", rect_b))

            if pattern == 0: # Clustered --------------------- no diffusion break ---------------------------------------
                # drain metal and via array
                for pnt in points:
                    loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
                    n0 = self.total_fingers // 2 + int(self.total_fingers % 2)

                    for i in range(con[0]):
                        add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'],
                                               pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext),
                                      self.m0.stack * 2 * lp['M2']['pitch'], nd0, 'H', "V1")

                        xi, yi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext, pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln0 = self.m0.stack*(2*(nd0-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d0 = [[xi, yi], [xi+ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR("d0", rect_d0))

                        ln1 = self.m1.stack*(2*(nd1-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d1 = [[xi+ln1-lp["M1"]["width"]-2*C.m_m_ext+2*self.m1.stack*lp['M1']['pitch'], yi],[xi+ln1-lp["M1"]["width"]-2*C.m_m_ext+2*lp['M1']['pitch'] +ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        add_via_array(cell_i, (rect_d1[0][0]+C.m_m_ext, rect_d1[0][1]), self.m0.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")
                        contact_rects.append(CR("d1", rect_d1))

            elif pattern == 2:  # Interdigitated -------------- no diffusion break ----------------------------------------
                for pnt in points:
                    nd0 = self.m0.fingers//2 + self.m0.fingers % 2
                    nd1 = self.m1.fingers//2 + self.m0.fingers % 2

                    for i in range(con[0]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        yi = pnt["p0"][0][1] + v_sp + (con[1]+i) * lp["M2"]['pitch'] + C.m_m_ext
                        ln0 = self.m0.stack*(4*(nd0-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d0 = [[xi, yi], [xi+ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR("d0", rect_d0))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m0.stack * 4 * lp['M2']['pitch'], nd0, 'H', "V1")

                    for i in range(con[1]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext + 2*lp["M2"]["pitch"]
                        yi = pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln1 = self.m1.stack*(4*(nd1-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d1 = [[xi, yi], [xi+ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        contact_rects.append(CR("d1", rect_d1))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m1.stack * 4 * lp['M2']['pitch'], nd1, 'H', "V1")

            elif pattern == 4:  # Common Centroid ------------- no diffusion break --------------------------------------
                lbl = np.array(["d0", "d1"])
                for k, pnt in enumerate(points):
                    nd0 = self.m0.fingers//2 + self.m0.fingers % 2
                    nd1 = self.m1.fingers//2 + self.m0.fingers % 2
                    lbl = lbl[::-1]

                    for i in range(con[0]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        yi = pnt["p0"][0][1] + v_sp + (con[1]+i) * lp["M2"]['pitch'] + C.m_m_ext
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'], nd0//2, 'H', "V1")
                        x_gap = self.stack*(self.m0.fingers//2 + self.m1.fingers)*lp["M1"]["pitch"]
                        add_via_array(cell_i, (xi+x_gap+C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'], nd0//2, 'H', "V1")

                        ln0 = rect_s[1][0] - rect_s[0][0]
                        rect_d0 = [[rect_s[0][0], yi], [rect_s[0][0]+ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR(lbl[0], rect_d0))

                    x_init = self.stack * (self.m0.fingers // 2) * lp["M1"]["pitch"]
                    xi = xi + x_init
                    for i in range(con[1]):
                        yi = pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln1 = rect_s[1][0] - rect_s[0][0]
                        rect_d1 = [[rect_s[0][0], yi], [rect_s[0][0]+ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        contact_rects.append(CR(lbl[1], rect_d1))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m1.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")



        elif pattern == 1 or pattern == 3 or pattern == 5:      # diffusion break
            contact_rects = []

            if pattern == 1:   # Clustered -------------------------------------------------
                cell_0, c_p0 = create_mos(self.m0, con=[con[0], 1, con[2], 1], orientation="V", fabric_on=False)
                cell_1, c_p1 = create_mos(self.m1, con=[con[1], 1, con[2], 1], orientation="V", fabric_on=False)

                x_offset = (2*lp['poly']['dummies'] + self.stack*self.m0.fingers - 1) * lp['poly']['pitch']
                cell_i = create_empty_cell("cm", 1e-9, 1e-12)
                add_transformed_polygons(cell_0, cell_i, (0, 0))
                add_transformed_polygons(cell_1, cell_i, (x_offset, 0))

                for c in c_p0:
                    c.translate(0, 0)
                    if c.id == "d":
                        c.id = "d0"
                    contact_rects.append(c)
                for c in c_p1:
                    c.translate(x_offset, 0)
                    if c.id == "d":
                        c.id = "d1"
                    contact_rects.append(c)

                for i in range(len(contact_rects)):
                    for j in range(i+1, len(contact_rects)):
                        ci = contact_rects[i]
                        cj = contact_rects[j]

                        if ((ci.rect[1][0] - ci.rect[0][0]) > (ci.rect[1][1] - ci.rect[0][1])) and ((cj.rect[1][0] - cj.rect[0][0]) > (cj.rect[1][1] - cj.rect[0][1])):
                            if set([ci.rect[0][1], ci.rect[1][1]]) == set([cj.rect[0][1], cj.rect[1][1]]) and ci.rect[1][0]< cj.rect[0][0]:
                                bridge_rect = [[ci.rect[1][0]-10, ci.rect[0][1]], [cj.rect[0][0]+10, cj.rect[1][1]]]
                                add_metal(cell_i, P(bridge_rect[0][0],bridge_rect[0][1]), bridge_rect[1][0]-bridge_rect[0][0], "H", "M2")
                                contact_rects.append(CR(ci.id, bridge_rect))



            elif pattern == 3:  # Interigitated ------------------------------------------
                m = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack, 'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x = create_empty_cell("temp", unit=1e-9, precision=1e-12)
                    points = create_base_layout(cell_x, m.fingers, m.fins, m.mos_type, stack=m.stack, multiplier=m.multiplier, orientation="V")

                    x_offset = (2 * lp['poly']['dummies'] + self.stack - 1) * lp['poly']['pitch']
                    add_transformed_polygons(cell_x, cell_i, (i * x_offset, 0))
                    lab_ = ["d0", "d1"]
                    con_ = [con[0], con[1]]

                    for p, pnt in enumerate(points):
                        l = divide_line(pnt["p0"][0][1], pnt["p1"][0][1], con[0] + con[1] + con[2], C.m_m_ext, "M2")
                        ln = (self.stack)*lp["M2"]["pitch"] + lp["M2"]["width"]
                        xi = pnt["s"][0][0]-C.m_m_ext + i*x_offset

                        for j in range(con_[p % 2]):
                            rect_d0 = [[xi, l[j]], [xi + ln+2*C.m_m_ext, l[j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[j]), ln+2*C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[p % 2], rect_d0))
                            if i % 2:
                                add_via(cell_i, [pnt["d"][0][0] + i*x_offset, l[j]], "V1")

                        for j in range(con_[(p+1) % 2]):
                            rect_d1 = [[xi, l[con_[p % 2] + j]], [xi + ln+2*C.m_m_ext, l[con_[p % 2] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con_[p % 2] + j]), ln+2*C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[(p+1) % 2], rect_d1))
                            if (i+1) % 2:
                                add_via(cell_i, [pnt["d"][0][0] + i*x_offset, l[con[1] + j]], "V1")

                        for j in range(con[2]):
                            rect_ds = [[xi, l[con[1] + con[0] + j]], [xi + ln+2*C.m_m_ext, l[con[1] + con[0] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con[1] + con[0] + j]), ln+2*C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("s", rect_ds))
                            add_via(cell_i, [pnt["s"][0][0] + i*x_offset, l[con[1] + con[0] + j]], "V1")

                        # gate
                        xi = pnt["s"][0][0] - C.m_m_ext + i*x_offset
                        yi = pnt["g"][0][1]
                        rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                        add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                        contact_rects.append(CR("g", rect_g))
                        rect = [[pnt["g"][0][0] + i*x_offset, pnt["g"][0][1]], [pnt["g"][0][0] + self.stack*lp["M1"]["pitch"] - lp["M1"]["width"] + i*x_offset, pnt["g"][0][1] + 32]]
                        fill_area_vias(cell_i, rect, "V1")

                        # bulk
                        if pnt["b"]:
                            xi = pnt["s"][0][0] - C.m_m_ext + i*x_offset
                            yi = pnt["b"][0][1] + 4
                            rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("b", rect_g))
                            add_via(cell_i, [pnt["b"][0][0]+i*x_offset+24, yi], "V1")
                #
                for i in range(len(contact_rects)):
                    for j in range(i+1, len(contact_rects)):
                        ci = contact_rects[i]
                        cj = contact_rects[j]
                        if ((ci.rect[1][0] - ci.rect[0][0]) > (ci.rect[1][1] - ci.rect[0][1])) and ((cj.rect[1][0] - cj.rect[0][0]) > (cj.rect[1][1] - cj.rect[0][1])):
                            if set([ci.rect[0][1], ci.rect[1][1]]) == set([cj.rect[0][1], cj.rect[1][1]]) and ci.rect[1][0]< cj.rect[0][0]:
                                bridge_rect = [[ci.rect[1][0]-10, ci.rect[0][1]], [cj.rect[0][0]+10, cj.rect[1][1]]]
                                add_metal(cell_i, P(bridge_rect[0][0],bridge_rect[0][1]), bridge_rect[1][0]-bridge_rect[0][0], "H", "M2")
                                contact_rects[i].rect[1][0] = contact_rects[j].rect[1][0]
                                contact_rects[j].id = "x"

                                # contact_rects.pop(j)
                                contact_rects.append(CR(ci.id, bridge_rect))



            elif pattern == 5:  # Common-Centroid -----------------------------------------
                m0 = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack,
                          'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x, c_pi = create_mos(m0, con=[0, 0, 0, 0], orientation="V", fabric_on=False)
                    add_transformed_polygons(cell_x, cell_i, (i * (2 * lp['poly']['dummies'] + self.stack * 1 - 1) * lp['poly']['pitch'], 0))



        # for cr in contact_rects:
        #     print(cr.id, cr.rect)
        #
        # c_s, s_s   = [[0,0],[0,0]],[[0,0],[0,0]]
        # c_d0, s_d0 = [[0,0],[0,0]],[[0,0],[0,0]]
        # c_d1, s_d1 = [[0,0],[0,0]],[[0,0],[0,0]]
        # c_b, s_b   = [[0,0],[0,0]],[[0,0],[0,0]]
        # c_g, s_g   = [[0,0],[0,0]],[[0,0],[0,0]]
#
# # -------------- Adding Labels ---------------------------------------------
        rect_s = []
        rect_b = []
        rect_d0 = []
        rect_d1 = []
        rect_g = []
        for cr in contact_rects:
            if (cr.rect[1][0] - cr.rect[0][0]) > (cr.rect[1][1] - cr.rect[0][1]):  # check if it's horizontal
                if cr.id == "s":
                    rect_s.append(cr.rect)
                elif cr.id == "d0":
                    rect_d0.append(cr.rect)
                elif cr.id == "d1":
                    rect_d1.append(cr.rect)
                elif cr.id == "b":
                    rect_b.append(cr.rect)
                elif cr.id == "g":
                    rect_g.append(cr.rect)

        c_s,  s_s = find_center_n_span_rects(rect_s)
        c_d0, s_d0 = find_center_n_span_rects(rect_d0)
        c_d1, s_d1 = find_center_n_span_rects(rect_d1)
        c_b, s_b = find_center_n_span_rects(rect_b)
        c_g, s_g = find_center_n_span_rects(rect_g)
        if labels[2]:
            add_label(cell_i, labels[2], P(c_s[0], c_s[1]), "M2LABEL")
        if labels[1]:
            add_label(cell_i, labels[1], P(c_d1[0], c_d1[1]), "M2LABEL")
        if labels[0]:
            add_label(cell_i, labels[0], P(c_d0[0], c_d0[1]), "M2LABEL")
        #--------------------------------------------------------------------------------------

        #------- Adding the Vertical M3 lines -------------------------------------------------
        if pattern == 0:
            h_sp = (con[2]-1)*lp["M3"]["pitch"]
            for i in range(con[2]):
                xi = c_s[0] - h_sp/2 - lp["M3"]["width"]/2 + i*lp["M3"]["pitch"]
                rect = [(xi, s_s[0][1]-C.m_m_ext), (xi+lp["M3"]["width"], s_b[1][1]+C.m_m_ext)]
                add_metal(cell_i, P(xi, s_s[0][1]-C.m_m_ext), s_b[1][1] - s_s[0][1] + 2*C.m_m_ext, "V", "M3")
                contact_rects.append(CR("s", rect))

            h_sp = (con[1]-1)*lp["M3"]["pitch"]
            for i in range(con[1]):
                xi = c_d0[0] - h_sp/2 - lp["M3"]["width"]/2 + i*lp["M3"]["pitch"]
                rect = [(xi, s_d0[0][1]-C.m_m_ext), (xi+lp["M3"]["width"], s_g[1][1] + C.m_m_ext)]
                add_metal(cell_i, P(rect[0][0], rect[0][1]), rect[1][1]-rect[0][1], "V", "M3")
                contact_rects.append(CR("d0", rect))

            h_sp = (con[0] - 1) * lp["M3"]["pitch"]
            for i in range(con[0]):
                xi = c_d1[0] - h_sp / 2 - lp["M3"]["width"] / 2 + i * lp["M3"]["pitch"]
                rect = [(xi, s_d1[0][1] - C.m_m_ext), (xi+lp["M3"]["width"], s_d1[1][1] + C.m_m_ext)]
                add_metal(cell_i, P(xi, s_d1[0][1] - C.m_m_ext), s_d1[1][1] - s_d1[0][1] + 2 * C.m_m_ext, "V", "M3")
                contact_rects.append(CR("d1", rect))

            # adding vias at the intersections
            for i in range(len(contact_rects)):
                for j in range(i+1, len(contact_rects)):
                    if(contact_rects[i].id == contact_rects[j].id):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["d0", "g"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["s", "b"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

        if pattern == 2 or pattern == 4:
            h_sp = (con[0] + con[1] + con[2] - 1)*lp["M3"]["pitch"]

            for i in range(con[0]):
                xi = c_s[0] - h_sp / 2 - lp["M3"]["width"] / 2 + i * lp["M3"]["pitch"]
                rect = [(xi, s_d0[0][1] - C.m_m_ext), (xi + lp["M3"]["width"], s_g[1][1] + C.m_m_ext)]
                add_metal(cell_i, P(xi, s_d0[0][1] - C.m_m_ext), s_g[1][1] - s_d0[0][1] + 2 * C.m_m_ext, "V", "M3")
                contact_rects.append(CR("d0", rect))

            x_i = xi + lp["M3"]["pitch"]
            for i in range(con[2]):
                xi = x_i + i*lp["M3"]["pitch"]
                rect = [(xi, s_s[0][1] - C.m_m_ext), (xi + lp["M3"]["width"], s_b[1][1] + C.m_m_ext)]
                add_metal(cell_i, P(xi, s_s[0][1] - C.m_m_ext), s_b[1][1] - s_s[0][1] + 2 * C.m_m_ext, "V", "M3")
                contact_rects.append(CR("s", rect))

            x_i += con[2]*lp["M3"]["pitch"]
            for i in range(con[1]):
                xi = x_i + i*lp["M3"]["pitch"]
                rect = [(xi, s_d1[0][1] - C.m_m_ext), (xi + lp["M3"]["width"], s_d1[1][1] + C.m_m_ext)]
                add_metal(cell_i, P(rect[0][0], rect[0][1]), rect[1][1] - rect[0][1], "V", "M3")
                contact_rects.append(CR("d1", rect))


            # adding vias at the intersections
            for i in range(len(contact_rects)):
                for j in range(i + 1, len(contact_rects)):
                    if (contact_rects[i].id == contact_rects[j].id):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["d0", "g"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["s", "b"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

        if pattern == 1:
            for cr in contact_rects:
                rect_d_g = [[0,0],[0,0]]
                rect_s_b = [[0,0],[0,0]]
                if cr.id == "d0" and cr.metal == "M3":
                    rect_d_g = [[cr.rect[0][0], cr.rect[1][1] - C.mm_overlap], [cr.rect[0][0]+lp["M3"]["width"], s_g[1][1] + C.mm_overlap + C.m_m_ext]]
                    add_metal(cell_i, P(rect_d_g[0][0], rect_d_g[0][1]), rect_d_g[1][1]-rect_d_g[0][1], "V", "M3")

                elif cr.id == "s" and cr.metal == "M3":
                    rect_s_b = [[cr.rect[0][0], cr.rect[1][1] - C.mm_overlap], [cr.rect[0][0]+lp["M3"]["width"], s_b[1][1] + C.mm_overlap + C.m_m_ext]]
                    add_metal(cell_i, P(rect_s_b[0][0], rect_s_b[0][1]), rect_s_b[1][1]-rect_s_b[0][1], "V", "M3")

                for cc in contact_rects:
                     if cc.metal == "M2" and cc.id == "g":
                         add_vias_at_recti_v_rectj(cell_i, rect_d_g, cc.rect, "M2")
                     elif cc.metal == "M2" and cc.id == "b":
                         add_vias_at_recti_v_rectj(cell_i, rect_s_b, cc.rect, "M2")

        if pattern == 3:
            p = cell_i.bounding_box()
            x_l = divide_line(p[0][0], p[1][0], con[0]+con[1]+con[2], 20, "M3")

            for i in range(con[0]):
                rd0 = [[x_l[i], s_d0[0][1]-C.m_m_ext], [x_l[i] + lp["M2"]["width"], s_g[1][1]+C.m_m_ext]]
                add_metal(cell_i, P(rd0[0][0], rd0[0][1]), rd0[1][1] - rd0[0][1], "V", "M3")
                contact_rects.append(CR("d0", rd0, "M3"))

            for i in range(con[2]):
                rs = [[x_l[con[0]+i], s_s[0][1]-C.m_m_ext], [x_l[con[0]+i] + lp["M2"]["width"], s_b[1][1]+C.m_m_ext]]
                add_metal(cell_i, P(rs[0][0], rs[0][1]), rs[1][1] - rs[0][1], "V", "M3")
                contact_rects.append(CR("s", rs, "M3"))

            for i in range(con[1]):
                rd1 = [[x_l[con[0]+con[2]+i], s_d1[0][1]-C.m_m_ext], [x_l[con[0]+con[2]+i]+lp["M2"]["width"], s_d1[1][1]+C.m_m_ext]]
                print("yes:", s_d1[1][1]-s_d1[0][1])
                add_metal(cell_i, P(rd1[0][0], rd1[0][1]), rd1[1][1]-rd1[0][1], "V", "M3")
                contact_rects.append(CR("d1", rd1, "M3"))


            # adding vias at the intersections
            for i in range(len(contact_rects)):
                for j in range(i + 1, len(contact_rects)):
                    if ((contact_rects[i].id == contact_rects[j].id) and
                            (contact_rects[i].metal == "M3" or contact_rects[j].metal == "M3")):

                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["d0", "g"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

                    elif set([contact_rects[i].id, contact_rects[j].id]) == set(["s", "b"]):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")
        # -------------------------------------------------------------------------------------
















        # -------- poly and fins fabric inclusion ------------------------------------------
        p = cell_i.bounding_box()
        poly_n = math.ceil((p[1][0] - p[0][0]) / lp["poly"]["pitch"]) - 1
        fin_n = math.ceil((p[1][1] - p[0][1]) // lp["fins"]["pitch"])
        x_ = 0
        y_ = 0
        if fabric_on:
            create_fabric(self.cell, poly_n, fin_n)
            x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
            y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width']) / 2

        add_transformed_polygons(cell_i, self.cell, (x_, y_))
