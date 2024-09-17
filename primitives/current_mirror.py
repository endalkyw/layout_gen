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

    def pin_polygon(self, rect: list, label: str, metal: str):
        pass

    def add_labels(self):
        pass

    def create_layout(self, pattern: int,
                      labels = ["do", "d1", "s"], label_pos = ["c","c","c"],
                      con = [1, 1, 1], fabric_on=True):

        assert pattern == 0 or 1 or 2 or 3 or 4 or 5
        fin_no = self.total_fins + 4 * lp['fins']['dummies'] - 2
        poly_no = self.stack*self.total_fingers + 2 * lp['poly']['dummies']
        cell_i = create_empty_cell("cm", 1e-9, 1e-12)

        n0 = self.m0.fingers//2 + int(self.m0.fingers % 2)
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
                        contact_rects.append(CR("d", rect_d0))

                        ln1 = self.m1.stack*(2*(nd1-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d1 = [[xi+ln1-lp["M1"]["width"]-2*C.m_m_ext+2*self.m1.stack*lp['M1']['pitch'], yi],[xi+ln1-lp["M1"]["width"]-2*C.m_m_ext+2*lp['M1']['pitch'] +ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        add_via_array(cell_i, (rect_d1[0][0]+C.m_m_ext, rect_d1[0][1]), self.m0.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")
                        contact_rects.append(CR("d", rect_d0))

            elif pattern == 2:  # Interdigitated -------------- no diffusion break ----------------------------------------
                for pnt in points:
                    loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
                    n0 = self.total_fingers // 2 + int(self.total_fingers % 2)

                    nd0 = self.m0.fingers//2 + self.m0.fingers % 2
                    nd1 = self.m1.fingers//2 + self.m0.fingers % 2

                    for i in range(con[0]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        yi = pnt["p0"][0][1] + v_sp + (con[1]+i) * lp["M2"]['pitch'] + C.m_m_ext
                        ln0 = self.m0.stack*(4*(nd0-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d0 = [[xi, yi], [xi+ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR("d", rect_d0))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m0.stack * 4 * lp['M2']['pitch'], nd0, 'H', "V1")

                    for i in range(con[1]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext + 2*lp["M2"]["pitch"]
                        yi = pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln1 = self.m1.stack*(4*(nd1-1))*lp["M1"]["pitch"] + lp["M1"]["width"] + 2*C.m_m_ext
                        rect_d1 = [[xi, yi], [xi+ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        contact_rects.append(CR("d", rect_d0))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m1.stack * 4 * lp['M2']['pitch'], nd1, 'H', "V1")

            elif pattern == 4:  # Common Centroid ------------- no diffusion break --------------------------------------
                for pnt in points:
                    loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
                    n0 = self.total_fingers // 2 + int(self.total_fingers % 2)

                    nd0 = self.m0.fingers//2 + self.m0.fingers % 2
                    nd1 = self.m1.fingers//2 + self.m0.fingers % 2

                    for i in range(con[0]):
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        yi = pnt["p0"][0][1] + v_sp + (con[1]+i) * lp["M2"]['pitch'] + C.m_m_ext
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'], nd0//2, 'H', "V1")
                        x_gap = self.stack*(self.m0.fingers//2 + self.m1.fingers)*lp["M1"]["pitch"]
                        add_via_array(cell_i, (xi+x_gap+C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'], nd0//2, 'H', "V1")

                        ln0 = rect_s[1][0] - rect_s[0][0]
                        rect_d0 = [[rect_s[0][0], yi], [rect_s[0][0]+ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR("d", rect_d0))
                    x_init = self.stack * (self.m0.fingers // 2) * lp["M1"]["pitch"]
                    xi = xi + x_init
                    for i in range(con[1]):
                        yi = pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln1 = rect_s[1][0] - rect_s[0][0]
                        rect_d1 = [[rect_s[0][0], yi], [rect_s[0][0]+ln1, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        contact_rects.append(CR("d", rect_d0))
                        add_via_array(cell_i, (xi+C.m_m_ext, yi), self.m1.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")




        elif pattern == 1 or pattern == 3 or pattern == 5:      # diffusion break
            contact_rects = []

            if pattern == 1:   # Clustered -------------------------------------------------
                cell_0, c_p0 = create_mos(self.m0, con=[con[0], 1, con[2], 1], orientation="V", fabric_on=False)
                cell_1, c_p1 = create_mos(self.m1, con=[con[1], 1, con[2], 1], orientation="V", fabric_on=False)

                cell_i = create_empty_cell("cm", 1e-9, 1e-12)
                add_transformed_polygons(cell_0, cell_i, (0, 0))
                add_transformed_polygons(cell_1, cell_i, ((2*lp['poly']['dummies'] + self.stack*self.m0.fingers - 1) * lp['poly']['pitch'], 0))

            elif pattern == 3:  # Interdigitated ------------------------------------------
                m0 = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack, 'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x, c_pi = create_mos(m0, con=[con[0], 1, con[2], 1], orientation="V", fabric_on=False)
                    add_transformed_polygons(cell_x, cell_i, (i*(2*lp['poly']['dummies'] + self.stack*1 - 1) * lp['poly']['pitch'], 0))

            elif pattern == 5:  # Common-Centroid -----------------------------------------
                m0 = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack,
                          'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x, c_pi = create_mos(m0, con=[0, 0, 0, 0], orientation="V", fabric_on=False)
                    add_transformed_polygons(cell_x, cell_i, (
                    i * (2 * lp['poly']['dummies'] + self.stack * 1 - 1) * lp['poly']['pitch'], 0))









        # poly and fins fabric inclusion ------------------------------------------
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
