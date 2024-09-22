from utilities.elements_utils import *
from core import create_base_layout, create_fabric
from utilities.pvector import P
from utilities.via import *
from primitives.mos import *

from core.constants import Const as C
class differential_pair:
    def __init__(self, m0: Mos, m1: Mos, cell_name):
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
        self.s = P(0, 0)
        self.mos_type = m0.mos_type
        self.pin_positions = {"d0": [], "d1": [], "s": []}

        self.lib = gdstk.Library(name=cell_name, unit=1e-9, precision=1e-12)
        self.cell = self.lib.new_cell(cell_name)

    def create_layout(self, pattern: int,
                      labels=["b", "g0", "g1", "d0", "d1", "s"], label_pos=["c", "c", "c"],
                      con=[1, 1, 1], fabric_on=True):

        assert pattern == 0 or 1 or 2 or 3 or 4 or 5
        cell_i = create_empty_cell("dp", 1e-9, 1e-12)

        nd0 = self.m0.fingers // 2 + self.m0.fingers % 2
        nd1 = self.m1.fingers // 2 + self.m0.fingers % 2

        if pattern == 0:
            g_t = "dp_0"
        elif pattern == 2:
            g_t = "dp_1"
        elif pattern == 4:
            g_t = "dp_2"

        if pattern == 0 or pattern == 2 or pattern == 4:  # no diffusion break
            # common to all connections: source, gate, and the bulk is common (no need to add connection)
            points = create_base_layout(cell_i, self.total_fingers, self.total_fins, self.m0.mos_type,
                                        stack=self.m0.stack, multiplier=self.m0.multiplier, orientation="V",
                                        body_contact=False, gate_type=g_t)
            contact_rects = []
            print(points)


            for pnt in points:
                loc = {"x_span": [pnt["p0"][0][0], pnt["p1"][0][0]], "y_span": [pnt["p0"][0][1], pnt["p1"][0][1]]}
                # metal connections
                if pattern == 0:
                    src_drn_cons = con[0] + con[2]
                    drn = con[0]

                else:
                    src_drn_cons = con[0] + con[1] + con[2]
                    drn = con[0] + con[1]

                v_pos = divide_line(loc["y_span"][0],loc["y_span"][1], src_drn_cons, C.m_m_ext, "M2")

                avail_v_space = loc["y_span"][1] - loc["y_span"][0]
                sd_space = (src_drn_cons - 1) * lp["M2"]["pitch"] + lp["M2"]["width"] + 2 * C.m_m_ext
                v_sp = (avail_v_space - sd_space) / 2

                # source metal and via array
                for i in range(con[2]):
                    xi, yi = loc["x_span"][0], pnt["p0"][0][1] + v_sp + drn * lp["M2"]['pitch'] + i * lp["M2"][
                        'pitch'] + C.m_m_ext
                    rect_s = [[xi, yi], [xi + loc["x_span"][1] - loc["x_span"][0], yi + lp["M2"]['width']]]
                    add_metal(cell_i, P(xi, yi), loc["x_span"][1] - loc["x_span"][0], "H", "M2")
                    add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1,
                                           pnt["p0"][0][1] + v_sp + drn * lp["M2"]['pitch'] + i * lp["M2"][
                                               'pitch'] + C.m_m_ext), self.m0.stack * 2 * lp['M1']['pitch'],
                                  self.total_fingers // 2 + 1, 'H', "V1")
                    contact_rects.append(CR("s", rect_s))

                # bulk metal and via array
                if pnt["b"]:
                    rect_b = [(loc["x_span"][0] + C.rx_m1, pnt["b"][0][1] + 4), (
                        loc["x_span"][0] + C.rx_m1 + loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1,
                        pnt["b"][0][1] + lp["M1"]["width"] + 4)]
                    add_vias_at_intersection(cell_i, rect_b, "M1")
                    add_metal(cell_i, P(rect_b[0][0], rect_b[0][1]), loc["x_span"][1] - loc["x_span"][0] - 2 * C.rx_m1,
                              "H",
                              "M2")
                    contact_rects.append(CR("b", rect_b))

            if pattern == 0:  # Clustered --------------------- no diffusion break ---------------------------------------
                # drain metal and via array
                for pnt in points:
                    # gate connection
                    g_n = ["g0", "g1"]
                    for i, g in enumerate(pnt["g"]):
                        rect_g = [[g[0]+C.m_m_ext, g[1]-1], [g[0]+self.stack*78*self.m0.fingers-32-C.m_m_ext, g[1]+32]]
                        add_vias_at_intersection(cell_i, rect_g, "M1")
                        add_metal(cell_i, P(rect_g[0][0], rect_g[0][1]),rect_g[1][0]-rect_g[0][0],"H","M2")
                        contact_rects.append(CR(g_n[i], rect_g))


                    for i in range(con[0]):
                        add_via_array(cell_i, (pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'],
                                               pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext),
                                      self.m0.stack * 2 * lp['M2']['pitch'], nd0, 'H', "V1")

                        xi, yi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext, \
                                 pnt["p0"][0][1] + v_sp + i * lp["M2"]['pitch'] + C.m_m_ext
                        ln0 = self.m0.stack * (2 * (nd0 - 1)) * lp["M1"]["pitch"] + lp["M1"]["width"] + 2 * C.m_m_ext
                        rect_d0 = [[xi, yi], [xi + ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR("d0", rect_d0))

                        ln1 = self.m1.stack * (2 * (nd1 - 1)) * lp["M1"]["pitch"] + lp["M1"]["width"] + 2 * C.m_m_ext
                        rect_d1 = [
                            [xi + ln1 - lp["M1"]["width"] - 2 * C.m_m_ext + 2 * self.m1.stack * lp['M1']['pitch'], yi],
                            [xi + ln1 - lp["M1"]["width"] - 2 * C.m_m_ext + 2 * self.m1.stack * lp['M1']['pitch'] + ln1,
                             yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln1, "H", "M2")
                        add_via_array(cell_i, (rect_d1[0][0] + C.m_m_ext, rect_d1[0][1]),
                                      self.m0.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")
                        contact_rects.append(CR("d1", rect_d1))

            elif pattern == 2:  # Interdigitated -------------- no diffusion break ----------------------------------------
                g_lst = [["g0", "g1"], ["g1", "g0"]]
                d_lst = [["d0", "d1"], ["d1", "d0"]]
                for ix, pnt in enumerate(points):
                    # gate connection
                    g_n = g_lst[ix % 2]
                    d_n = d_lst[ix % 2]

                    print(g_n[0], d_n[1])
                    for i, g in enumerate(sorted(pnt["g"])):
                        rect_g = [[min(pnt["g"][0][0], pnt["g"][1][0]) + C.m_m_ext, g[1] - 1],
                                  [min(pnt["g"][0][0],
                                       pnt["g"][1][0]) + self.stack * 78 * self.total_fingers - 32 - C.m_m_ext,
                                   g[1] + 32]]
                        add_vias_at_intersection(cell_i, rect_g, "M1")
                        add_metal(cell_i, P(rect_g[0][0], rect_g[0][1]), rect_g[1][0] - rect_g[0][0], "H", "M2")
                        contact_rects.append(CR(g_n[i], rect_g))

                    sp_v = divide_line(pnt["p0"][0][1], pnt["p1"][0][1], con[0] + con[1] + con[2], C.m_m_ext, "M2")
                    m1s = m1_array(pnt["p0"][0][0] + C.rx_m1, pnt["p1"][0][0], self.stack, 1)

                    i = 0
                    for _ in range(con[0]):  # d0
                        yi = sp_v[i]
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        ln0 = self.m0.stack * (self.total_fingers - 2) * lp["M1"]["pitch"] + lp["M1"][
                            "width"] + 2 * C.m_m_ext
                        rect_d0 = [[xi, yi], [xi + ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR(d_n[0], rect_d0))
                        add_via_array(cell_i, (m1s[1], yi), self.m0.stack * 4 * lp['M2']['pitch'], nd0,
                                      'H', "V1")
                        i += 1

                    for _ in range(con[1]):  # d1
                        yi = sp_v[i]
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        ln0 = self.m0.stack * (self.total_fingers - 2) * lp["M1"]["pitch"] + lp["M1"]["width"] + 2 * C.m_m_ext
                        rect_d1 = [[xi, yi], [xi + ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln0, "H", "M2")
                        add_via_array(cell_i, (m1s[3], yi), self.m1.stack * 4 * lp['M2']['pitch'], nd1, 'H', "V1")
                        contact_rects.append(CR(d_n[1], rect_d1))
                        i += 1

            elif pattern == 4:  # Common Centroid ------------- no diffusion break --------------------------------------
                g_lst = [["g0", "g1", "g0"], ["g1", "g0", "g1"]]
                d_lst = [["d0", "d1"], ["d1", "d0"]]
                for ix, pnt in enumerate(points):
                    # gate connection

                    g_n = g_lst[ix%2]
                    d_n = d_lst[ix%2]

                    for i, g in enumerate(sorted(pnt["g"])):
                        rect_g = [[min(pnt["g"][0][0], pnt["g"][1][0]) + C.m_m_ext, g[1] - 1],
                                  [min(pnt["g"][0][0],
                                       pnt["g"][1][0]) + self.stack * 78 * self.total_fingers - 32 - C.m_m_ext,
                                   g[1] + 32]]
                        add_vias_at_intersection(cell_i, rect_g, "M1")
                        add_metal(cell_i, P(rect_g[0][0], rect_g[0][1]), rect_g[1][0] - rect_g[0][0], "H", "M2")
                        contact_rects.append(CR(g_n[i], rect_g))

                    sp_v = divide_line(pnt["p0"][0][1], pnt["p1"][0][1], con[0]+con[1]+con[2], C.m_m_ext, "M2")
                    m1s  = m1_array(pnt["p0"][0][0]+C.rx_m1, pnt["p1"][0][0], self.stack, 1)

                    i = 0
                    for _ in range(con[0]): # d0
                        yi = sp_v[i]
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        ln0 = self.m0.stack * (self.total_fingers-2) * lp["M1"]["pitch"] + lp["M1"]["width"] + 2 * C.m_m_ext
                        rect_d0 = [[xi, yi], [xi + ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d0[0][0], rect_d0[0][1]), ln0, "H", "M2")
                        contact_rects.append(CR(d_n[0], rect_d0))

                        add_via_array(cell_i, (xi + C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'], nd0 // 2,'H', "V1")
                        x_gap = self.stack * (self.m0.fingers // 2 + self.m1.fingers) * lp["M1"]["pitch"]
                        add_via_array(cell_i, (xi + x_gap + C.m_m_ext, yi), self.m0.stack * 2 * lp['M2']['pitch'],nd0 // 2, 'H', "V1")
                        i += 1

                    for _ in range(con[1]): # d1
                        yi = sp_v[i]
                        xi = pnt["p0"][0][0] + C.rx_m1 + self.m0.stack * lp['M1']['pitch'] - C.m_v_ext
                        ln0 = self.m0.stack * (self.total_fingers - 2) * lp["M1"]["pitch"] + lp["M1"]["width"] + 2 * C.m_m_ext
                        rect_d1 = [[xi, yi], [xi + ln0, yi + lp["M1"]["width"]]]
                        add_metal(cell_i, P(rect_d1[0][0], rect_d1[0][1]), ln0, "H", "M2")
                        add_via_array(cell_i, (m1s[self.m0.fingers//2+1], yi), self.m1.stack * 2 * lp['M2']['pitch'], nd1, 'H', "V1")
                        contact_rects.append(CR(d_n[1], rect_d1))
                        i += 1



        elif pattern == 1 or pattern == 3 or pattern == 5:  # diffusion break
            contact_rects = []

            if pattern == 1:  # Clustered -------------------------------------------------
                cell_0, c_p0 = create_mos(self.m0, con=[con[0], 1, con[2], 1], orientation="V", fabric_on=False)
                cell_1, c_p1 = create_mos(self.m1, con=[con[1], 1, con[2], 1], orientation="V", fabric_on=False)

                x_offset = (2 * lp['poly']['dummies'] + self.stack * self.m0.fingers - 1) * lp['poly']['pitch']
                cell_i = create_empty_cell("cm", 1e-9, 1e-12)
                add_transformed_polygons(cell_0, cell_i, (0, 0))
                add_transformed_polygons(cell_1, cell_i, (x_offset, 0))

                for c in c_p0:
                    c.translate(0, 0)
                    if c.id == "d":
                        c.id = "d0"
                    if c.id == "g":
                        c.id = "g0"
                    contact_rects.append(c)

                for c in c_p1:
                    c.translate(x_offset, 0)
                    if c.id == "d":
                        c.id = "d1"
                    if c.id == "g":
                        c.id = "g1"
                    contact_rects.append(c)

                for i in range(len(contact_rects)):
                    for j in range(i + 1, len(contact_rects)):
                        ci = contact_rects[i]
                        cj = contact_rects[j]
                        if ci.id != "g0" and ci.id != "g1":
                            if ((ci.rect[1][0] - ci.rect[0][0]) > (ci.rect[1][1] - ci.rect[0][1])) and (
                                    (cj.rect[1][0] - cj.rect[0][0]) > (cj.rect[1][1] - cj.rect[0][1])):
                                if set([ci.rect[0][1], ci.rect[1][1]]) == set([cj.rect[0][1], cj.rect[1][1]]) and \
                                        ci.rect[1][0] < cj.rect[0][0]:
                                    bridge_rect = [[ci.rect[1][0] - 10, ci.rect[0][1]], [cj.rect[0][0] + 10, cj.rect[1][1]]]
                                    add_metal(cell_i, P(bridge_rect[0][0], bridge_rect[0][1]),
                                              bridge_rect[1][0] - bridge_rect[0][0], "H", "M2")
                                    contact_rects.append(CR(ci.id, bridge_rect))

            elif pattern == 3:  # Interigitated ------------------------------------------
                m = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack,
                         'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x = create_empty_cell("temp", unit=1e-9, precision=1e-12)
                    points = create_base_layout(cell_x, m.fingers, m.fins, m.mos_type, stack=m.stack,
                                                multiplier=m.multiplier, orientation="V")

                    x_offset = (2 * lp['poly']['dummies'] + self.stack - 1) * lp['poly']['pitch']
                    add_transformed_polygons(cell_x, cell_i, (i * x_offset, 0))
                    lab_ = ["d0", "d1"]
                    lab_g = ["g0", "g1"]
                    con_ = [con[0], con[1]]

                    for p, pnt in enumerate(points):
                        l = divide_line(pnt["p0"][0][1], pnt["p1"][0][1], con[0] + con[1] + con[2], C.m_m_ext, "M2")
                        ln = (self.stack) * lp["M2"]["pitch"] + lp["M2"]["width"]
                        xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset

                        for j in range(con_[p % 2]):
                            rect_d0 = [[xi, l[j]], [xi + ln + 2 * C.m_m_ext, l[j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[p % 2], rect_d0))
                            if i % 2:
                                add_via(cell_i, [pnt["d"][0][0] + i * x_offset, l[j]], "V1")

                        for j in range(con_[(p + 1) % 2]):
                            rect_d1 = [[xi, l[con_[p % 2] + j]],
                                       [xi + ln + 2 * C.m_m_ext, l[con_[p % 2] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con_[p % 2] + j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[(p + 1) % 2], rect_d1))
                            if (i + 1) % 2:
                                add_via(cell_i, [pnt["d"][0][0] + i * x_offset, l[con[1] + j]], "V1")

                        for j in range(con[2]):
                            rect_ds = [[xi, l[con[1] + con[0] + j]],
                                       [xi + ln + 2 * C.m_m_ext, l[con[1] + con[0] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con[1] + con[0] + j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("s", rect_ds))
                            add_via(cell_i, [pnt["s"][0][0] + i * x_offset, l[con[1] + con[0] + j]], "V1")

                        # gate
                        xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset
                        yi = pnt["g"][0][1]
                        rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                        add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                        contact_rects.append(CR("g", rect_g))
                        rect = [[pnt["g"][0][0] + i * x_offset, pnt["g"][0][1]],
                                [pnt["g"][0][0] + self.stack * lp["M1"]["pitch"] - lp["M1"]["width"] + i * x_offset,
                                 pnt["g"][0][1] + 32]]
                        fill_area_vias(cell_i, rect, "V1")

                        # bulk
                        if pnt["b"]:
                            xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset
                            yi = pnt["b"][0][1] + 4
                            rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("b", rect_g))
                            add_via(cell_i, [pnt["b"][0][0] + i * x_offset + 24, yi], "V1")
                #
                for i in range(len(contact_rects)):
                    for j in range(i + 1, len(contact_rects)):
                        ci = contact_rects[i]
                        cj = contact_rects[j]
                        if ((ci.rect[1][0] - ci.rect[0][0]) > (ci.rect[1][1] - ci.rect[0][1])) and (
                                (cj.rect[1][0] - cj.rect[0][0]) > (cj.rect[1][1] - cj.rect[0][1])):
                            if set([ci.rect[0][1], ci.rect[1][1]]) == set([cj.rect[0][1], cj.rect[1][1]]) and \
                                    ci.rect[1][0] < cj.rect[0][0]:
                                bridge_rect = [[ci.rect[1][0] - 10, ci.rect[0][1]], [cj.rect[0][0] + 10, cj.rect[1][1]]]
                                add_metal(cell_i, P(bridge_rect[0][0], bridge_rect[0][1]),
                                          bridge_rect[1][0] - bridge_rect[0][0], "H", "M2")
                                contact_rects[i].rect[1][0] = contact_rects[j].rect[1][0]
                                contact_rects[j].id = "x"

                                # contact_rects.pop(j)
                                contact_rects.append(CR(ci.id, bridge_rect))

            elif pattern == 5:  # Common-Centroid -----------------------------------------
                m = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'stack': self.m0.stack,
                         'multiplier': self.m0.multiplier, 'mos_type': self.m0.mos_type})

                for i in range(self.total_fingers):
                    cell_x = create_empty_cell("temp", unit=1e-9, precision=1e-12)
                    points = create_base_layout(cell_x, m.fingers, m.fins, m.mos_type, stack=m.stack,
                                                multiplier=m.multiplier, orientation="V")

                    x_offset = (2 * lp['poly']['dummies'] + self.stack - 1) * lp['poly']['pitch']
                    add_transformed_polygons(cell_x, cell_i, (i * x_offset, 0))
                    lab_ = ["d0", "d1"]
                    con_ = [con[0], con[1]]

                    for p, pnt in enumerate(points):
                        l = divide_line(pnt["p0"][0][1], pnt["p1"][0][1], con[0] + con[1] + con[2], C.m_m_ext, "M2")
                        ln = (self.stack) * lp["M2"]["pitch"] + lp["M2"]["width"]
                        xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset

                        for j in range(con_[p % 2]):
                            rect_d0 = [[xi, l[j]], [xi + ln + 2 * C.m_m_ext, l[j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[p % 2], rect_d0))
                            if self.m0.fingers // 2 <= i < self.m0.fingers // 2 + self.m1.fingers:
                                add_via(cell_i, [pnt["d"][0][0] + i * x_offset, l[j]], "V1")

                        for j in range(con_[(p + 1) % 2]):
                            rect_d1 = [[xi, l[con_[p % 2] + j]],
                                       [xi + ln + 2 * C.m_m_ext, l[con_[p % 2] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con_[p % 2] + j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR(lab_[(p + 1) % 2], rect_d1))
                            if 0 <= i < self.m0.fingers // 2 or self.m1.fingers + self.m0.fingers // 2 <= i < self.total_fingers:
                                add_via(cell_i, [pnt["d"][0][0] + i * x_offset, l[con[1] + j]], "V1")

                        for j in range(con[2]):
                            rect_ds = [[xi, l[con[1] + con[0] + j]],
                                       [xi + ln + 2 * C.m_m_ext, l[con[1] + con[0] + j] + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, l[con[1] + con[0] + j]), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("s", rect_ds))
                            add_via(cell_i, [pnt["s"][0][0] + i * x_offset, l[con[1] + con[0] + j]], "V1")

                        # gate
                        xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset
                        yi = pnt["g"][0][1]
                        rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                        add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                        contact_rects.append(CR("g", rect_g))
                        rect = [[pnt["g"][0][0] + i * x_offset, pnt["g"][0][1]],
                                [pnt["g"][0][0] + self.stack * lp["M1"]["pitch"] - lp["M1"]["width"] + i * x_offset,
                                 pnt["g"][0][1] + 32]]
                        fill_area_vias(cell_i, rect, "V1")

                        # bulk
                        if pnt["b"]:
                            xi = pnt["s"][0][0] - C.m_m_ext + i * x_offset
                            yi = pnt["b"][0][1] + 4
                            rect_g = [[xi, yi], [xi + ln + 2 * C.m_m_ext, yi + lp["M2"]["width"]]]
                            add_metal(cell_i, P(xi, yi), ln + 2 * C.m_m_ext, "H", "M2")
                            contact_rects.append(CR("b", rect_g))
                            add_via(cell_i, [pnt["b"][0][0] + i * x_offset + 24, yi], "V1")

                for i in range(len(contact_rects)):
                    for j in range(i + 1, len(contact_rects)):
                        ci = contact_rects[i]
                        cj = contact_rects[j]
                        if ((ci.rect[1][0] - ci.rect[0][0]) > (ci.rect[1][1] - ci.rect[0][1])) and (
                                (cj.rect[1][0] - cj.rect[0][0]) > (cj.rect[1][1] - cj.rect[0][1])):
                            if set([ci.rect[0][1], ci.rect[1][1]]) == set([cj.rect[0][1], cj.rect[1][1]]) and \
                                    ci.rect[1][0] < cj.rect[0][0]:
                                bridge_rect = [[ci.rect[1][0] - 10, ci.rect[0][1]], [cj.rect[0][0] + 10, cj.rect[1][1]]]
                                add_metal(cell_i, P(bridge_rect[0][0], bridge_rect[0][1]),
                                          bridge_rect[1][0] - bridge_rect[0][0], "H", "M2")
                                contact_rects[i].rect[1][0] = contact_rects[j].rect[1][0]
                                contact_rects[j].id = "x"
                                contact_rects.append(CR(ci.id, bridge_rect))



        # -------------- Adding Labels ---------------------------------------------
        rect_s = []
        rect_b = []
        rect_d0 = []
        rect_d1 = []
        rect_g0 = []
        rect_g1 = []

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
                elif cr.id == "g0":
                    rect_g0.append(cr.rect)
                elif cr.id == "g1":
                    rect_g1.append(cr.rect)

        c_s, s_s = find_center_n_span_rects(rect_s)
        c_d0, s_d0 = find_center_n_span_rects(rect_d0)
        c_d1, s_d1 = find_center_n_span_rects(rect_d1)
        c_b, s_b = find_center_n_span_rects(rect_b)
        c_g0, s_g0 = find_center_n_span_rects(rect_g0)
        c_g1, s_g1 = find_center_n_span_rects(rect_g1)

        if labels[5]:
            add_label(cell_i, labels[5], P(c_s[0], c_s[1]), "M2LABEL")
        if labels[4]:
            add_label(cell_i, labels[4], P(c_d1[0], c_d1[1]), "M2LABEL")
        if labels[3]:
            add_label(cell_i, labels[3], P(c_d0[0], c_d0[1]), "M2LABEL")
        if labels[1]:
            add_label(cell_i, labels[1], P(c_g0[0], c_g0[1]), "M2LABEL")
        if labels[2]:
            add_label(cell_i, labels[2], P(c_g1[0], c_g1[1]), "M2LABEL")
        if labels[0]:
            add_label(cell_i, labels[0], P(c_b[0], c_b[1]), "M2LABEL")
        # --------------------------------------------------------------------------------------

        # ------- Adding the Vertical M3 lines -------------------------------------------------
        if pattern == 0:
            sp_left = divide_line(s_d0[0][0], s_d0[1][0], con[0]+1, C.m_m_ext, "M3")
            sp_center = divide_line(s_s[0][0], s_s[1][0], con[2], C.m_m_ext, "M3")
            sp_right = divide_line(s_d1[0][0], s_d1[1][0], con[1]+1, C.m_m_ext, "M3")
            for i, xi in enumerate(sp_left):
                if i == len(sp_left)//2:
                    r = add_metal(cell_i, P(xi, s_g0[0][1] - C.m_m_ext), s_g0[1][1] - s_g0[0][1] + 2 * C.m_m_ext, "V", "M3", True)
                    contact_rects.append(CR("g0", r))
                else:
                    r = add_metal(cell_i, P(xi, s_d0[0][1] - C.m_m_ext), s_d0[1][1] - s_d0[0][1] + 2 * C.m_m_ext, "V", "M3", True)
                    contact_rects.append(CR("d0", r))

            for i, xi in enumerate(sp_right):
                if i == len(sp_right)//2:
                    r = add_metal(cell_i, P(xi, s_g1[0][1] - C.m_m_ext), s_g1[1][1] - s_g1[0][1] + 2 * C.m_m_ext, "V", "M3", True)
                    contact_rects.append(CR("g1", r))
                else:
                    r = add_metal(cell_i, P(xi, s_d1[0][1] - C.m_m_ext), s_d1[1][1] - s_d1[0][1] + 2 * C.m_m_ext, "V", "M3", True)
                    contact_rects.append(CR("d1", r))

            for i, xi in enumerate(sp_center):
                r = add_metal(cell_i, P(xi, s_s[0][1] - C.m_m_ext), s_s[1][1] - s_s[0][1] + 2 * C.m_m_ext, "V", "M3", True)
                contact_rects.append(CR("s", r))

            # adding vias at the intersections
            for i in range(len(contact_rects)):
                for j in range(i + 1, len(contact_rects)):
                    if (contact_rects[i].id == contact_rects[j].id):
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

        if pattern == 4 or pattern == 2:
            x_sp = divide_line(s_d0[0][0],s_d0[1][0], con[0] + con[1] + con[2] + 2, C.m_m_ext, "M3")

            for i in range(con[0]): # d0
                r = add_metal(cell_i, P(x_sp[i], s_d0[0][1] - C.m_m_ext), s_d0[1][1] - s_d0[0][1] + 2 * C.m_m_ext, "V",
                              "M3", True)
                contact_rects.append(CR("d0", r, "M3"))

            i += 1
            r = add_metal(cell_i, P(x_sp[i], s_g1[0][1] - C.m_m_ext), s_g1[1][1] - s_g1[0][1] + 2 * C.m_m_ext, "V","M3", True)
            contact_rects.append(CR("g1", r, "M3"))


            for i in range(i+1, i+1+con[2]): # s
                r = add_metal(cell_i, P(x_sp[i], s_s[0][1] - C.m_m_ext), s_s[1][1] - s_s[0][1] + 2 * C.m_m_ext, "V",
                              "M3", True)
                contact_rects.append(CR("s", r, "M3"))

            i+=1
            r = add_metal(cell_i, P(x_sp[i], s_g0[0][1] - C.m_m_ext), s_g0[1][1] - s_g0[0][1] + 2 * C.m_m_ext, "V","M3", True)
            contact_rects.append(CR("g0", r, "M3"))

            for i in range(i+1, i+1+con[1]): # d1
                r = add_metal(cell_i, P(x_sp[i], s_d1[0][1] - C.m_m_ext), s_d1[1][1] - s_d1[0][1] + 2 * C.m_m_ext, "V",
                              "M3", True)
                contact_rects.append(CR("d1", r, "M3"))

            # adding vias at the intersections
            for i in range(len(contact_rects)):
                for j in range(i + 1, len(contact_rects)):
                    if (contact_rects[i].id == contact_rects[j].id) and "M3" in [contact_rects[i].metal, contact_rects[j].metal]:
                        add_vias_at_recti_v_rectj(cell_i, contact_rects[i].rect, contact_rects[j].rect, "M2")

        if pattern == 1:
            # only add the source vertical metal
            sp_s = divide_line(s_s[0][0], s_s[1][0], con[2], C.m_m_ext, "M3")
            for i in range(con[2]):  # d0
                r = add_metal(cell_i, P(sp_s[i], s_s[0][1] - C.m_m_ext), s_s[1][1] - s_s[0][1] + 2 * C.m_m_ext, "V",
                              "M3", True)
                contact_rects.append(CR("s", r, "M3"))
                for i in range(len(contact_rects)):
                    if(contact_rects[i].id == "s"):
                        add_vias_at_recti_v_rectj(cell_i, r, contact_rects[i].rect, "M2")


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
