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
                      labels = ["do", "d1", "s"],
                      label_pos = ["c","c","c"],
                      con = [1, 1, 1],
                      row_number = 1,
                      fabric_on=True):

        assert pattern == 0 or 1 or 2 or 3 or 4 or 5
        fin_no = self.total_fins + 4 * lp['fins']['dummies'] - 2
        poly_no = self.stack*self.total_fingers + 2 * lp['poly']['dummies']

        cell_i = create_empty_cell("cm", 1e-9, 1e-12)
        loc = create_base_layout(cell_i, self.total_fingers,
                                  self.total_fins,
                                  self.mos_type,
                                  stack=self.stack,
                                  gate_type="cm")

        n0 = self.m0.fingers//2 + int(self.m0.fingers % 2)
        n1 = self.m1.fingers//2 + int(self.m1.fingers % 2)
        k1 = C.rx_m1 + self.stack * lp["M1"]["pitch"]

        if pattern == 0 or pattern == 2 or pattern == 4:        # no diffusion break
            # common to all connections: source, gate, and the bulk is common (no need to add connection)
            y_s = 50
            y_g = self.total_fins*lp['fins']['pitch'] + 65

            # source metal and via array
            add_metal(cell_i, P(0, y_s), (self.stack*self.total_fingers+1)*lp['poly']['pitch'], "H", "M2")
            add_via_array(cell_i, (C.rx_m1, y_s), self.stack*2*lp['M1']['pitch'], self.total_fingers//2+1, 'H', "V1")
            self.pin_polygon([0, y_s, (self.stack*self.total_fingers+1)*lp['poly']['pitch']], "s", "M2")

            # gate metal
            add_metal(cell_i, P(0, y_g), (self.stack*self.total_fingers + 1) * lp['poly']['pitch'], "H", "M1")

            # source to bulk connection
            x_len = self.stack*self.total_fingers*lp["poly"]["pitch"]
            add_metal(cell_i, P(x_len/2, y_s), (self.total_fins+1)*lp["fins"]["pitch"]+208, "V", "M2")

            if pattern == 0: # Clustered --------------------- no diffusion break ---------------------------------------
                # drain 0 connection
                add_v_metal_array(cell_i, P(k1, -65), 66, self.stack*2*lp['M2']['pitch'], n0, "M1")
                add_via_array(cell_i, (k1, -60), self.stack*2*lp['M2']['pitch'], n0, 'H', "V1")
                add_metal(cell_i, P(k1 - 10, -60), self.stack*(n0-1)*2*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
                self.pin_polygon([k1 - 10, -60, self.stack*(n0-1)*2*lp['M2']['pitch']+lp['M1']['width']+20], "d0", "M2")

                # gate_drain
                add_v_metal_array(cell_i, P(k1, y_g-65), 66, self.stack*2*lp['M1']['pitch'], n0, "M1")

                # drain 1 connection
                k = self.stack*n0*2*lp['M2']['pitch']
                add_v_metal_array(cell_i, P(k1+k, -65), 66, self.stack*2*lp['M1']['pitch'], n1, "M1")
                add_via_array(cell_i, (k1+k, -60), self.stack*2*lp['M2']['pitch'], n1, 'H', "V1")
                add_metal(cell_i, P(k1+k - 10, -60), self.stack*(n1-1)*2*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
                self.pin_polygon([k1+k - 10, -60, self.stack*(n1-1)*2*lp['M2']['pitch']+lp['M1']['width']+20], "d1","M2")

            elif pattern == 2:  # Interdigitated -------------- no diffusion break ----------------------------------------
                # drain 0 connection (this can be avoided for optimization stage)
                add_v_metal_array(cell_i, P(k1, -65), 60,  self.stack*4*lp['M2']['pitch'], n0, "M1")
                add_via_array(cell_i, (k1, -60),  self.stack*4*lp['M2']['pitch'], n0, 'H', "V1")
                add_metal(cell_i, P(k1 - 10, -60),  self.stack*(n0-1)*4*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
                self.d0.update(k1+16, 16)
                # self.pin_polygon([k1 - 10, -60, self.stack*(n0-1)*2*lp['M2']['pitch']+lp['M1']['width']+20], "d0", "M2")

                # gate drain contact m0
                add_v_metal_array(cell_i, P(k1, y_g-66), 98, self.stack*4*lp['M1']['pitch'], n0, "M1")  # body drain contact m0

                # drain 1 connection
                k = self.stack*2*lp['M2']['pitch']
                add_v_metal_array(cell_i, P(k1+k, -65), 66, self.stack*4*lp['M1']['pitch'], n1, "M1")
                add_via_array(cell_i, (k1+k, -60), self.stack*4*lp['M2']['pitch'], n1, 'H', "V1")
                add_metal(cell_i, P(k1+k - 10, -60), self.stack*(n1-1)*4*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
                self.d1.update(k1+k + 16, -60 + 16)

            elif pattern == 4:  # Common Centroid ------------- no diffusion break --------------------------------------
                # drain 0_1 connection
                add_v_metal_array(cell_i, P(k1, -65), 66, self.stack*2*lp['M2']['pitch'], n0//2, "M1")
                add_via_array(cell_i, (k1, -66), self.stack*2*lp['M2']['pitch'], n0//2, 'H', "V1")
                self.d0.update(k1+16, -60+16)

                # gate drain contact m0
                add_v_metal_array(cell_i, P(k1, y_g-66), 98, self.stack*2*lp['M1']['pitch'], n0//2, "M1")  # body drain contact m0

                # drain 1__ connection
                k = self.stack*(n0/2)*2*lp['M2']['pitch']
                add_v_metal_array(cell_i, P(k1+k, -143), 143, self.stack*2*lp['M1']['pitch'], n1, "M1")
                add_via_array(cell_i, (k1+k, -143), self.stack*2*lp['M2']['pitch'], n1, 'H', "V1")
                add_metal(cell_i, P(k1+k - 10, -138), self.stack*(n1-1)*2*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
                self.d1.update(k1+k + 16, -138 + 16)

                # drain 0_2 connection
                k = self.stack*n0*3*lp['M2']['pitch']
                add_v_metal_array(cell_i, P(k1+k, -65), 66, self.stack*2*lp['M2']['pitch'], n0//2, "M1")
                add_via_array(cell_i, (k1+k, -60), self.stack*2*lp['M2']['pitch'], int(n0/2), 'H', "V1")

                # gate drain contact m0
                add_v_metal_array(cell_i, P(k1+k, y_g-66), 98, self.stack*2*lp['M1']['pitch'], n0//2, "M1")  # body drain contact m0
                add_metal(cell_i, P(k1 - 10, -60), self.stack*2*(n1+n0-1)*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")

            x_ = 0
            y_ = 0
            if fabric_on:
                create_fabric(self.cell, poly_no, fin_no)
                x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

            add_transformed_polygons(cell_i, self.cell, (x_, y_))

            # self.s.transform_p(x_, y_)
            # self.d0.transform_p(x_, y_)
            # self.d1.transform_p(x_, y_)
            #
            # add_label(self.cell, labels[0], self.d0, "M1LABEL")
            # add_label(self.cell, labels[1], self.d1, "M1LABEL")
            # add_label(self.cell, labels[2], self.s, "M1LABEL")

            self.add_labels()

        elif pattern == 1 or pattern == 3 or pattern == 5:      # diffusion break

            if pattern == 1: # Clustered -------------------------------------------------
                mos_a_cell = create_mos(self.m0, fabric_on = False)
                mos_b_cell = create_mos(self.m1, fabric_on = False)
                poly_no = self.total_fingers + 4 * lp['poly']['dummies']-1

                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                add_transformed_polygons(mos_a_cell, cell_i, (0, 0))
                add_transformed_polygons(mos_b_cell, cell_i, ((2*lp['poly']['dummies'] + self.m0.fingers - 1) * lp['poly']['pitch'] , 0))

                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions
                # 1) source connection
                add_metal(cell_i, P(self.m0.fingers*lp['poly']['pitch'], y_s), (2*lp['poly']['dummies'])*lp['poly']['pitch'], "H", "M2")
                self.s.update(23, y_s)
                self.s.transform_p(16, 16)

                # body drain contact m0  ** gate_drain?
                add_v_metal_array(cell_i, P(101, y_g-65), 97, 2*lp['M1']['pitch'], n0, "M1")

                # 2) drain is already connected

                # 3) gate connection
                add_metal(cell_i, P(self.m0.fingers*lp['poly']['pitch'], y_g), (2*lp['poly']['dummies'])*lp['poly']['pitch'], "H", "M1")

                # 4) source to bulk
                x_len = (self.m0.fingers + lp['poly']['dummies'])*lp["poly"]["pitch"]
                add_metal(cell_i, P(((self.m0.fingers+1)/2)*lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), (((self.total_fingers+2)/2)+4)*lp['poly']['pitch'], "H", "M2")
                add_metal(cell_i, P(x_len, y_s), (self.total_fins+1)*lp["fins"]["pitch"]+205, "V", "M2")

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))

                self.s.transform_p(x_, y_)
                add_label(self.cell, labels[2], self.s, "M1LABEL")

            elif pattern == 3: # Interdigitated ------------------------------------------
                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                poly_no = self.total_fingers*2*(lp['poly']['dummies'])+1
                for i in range(self.total_fingers):
                    mi = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'multiplier': 1, 'mos_type': self.m0.mos_type})
                    mos_i = create_mos(mi, fabric_on = False, con=[0,1,1,1])
                    add_transformed_polygons(mos_i, cell_i,(i*(2*lp['poly']['dummies']) * lp['poly']['pitch'] , 0))

                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions

                # 1) source connection
                add_metal(cell_i, P(0, y_s), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                self.s.update(23, y_s)
                self.s.transform_p(16, 16)

                # 2) drain connection
                k = 2*lp['poly']['dummies']*lp['poly']['pitch'] # gap between the mosfets
                add_v_metal_array(cell_i, P(101, -65), 60, 2*k, self.m0.fingers, "M1")
                add_via_array(cell_i, (101, -60), 2*k, self.m0.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16, -60), (self.m0.fingers-1)*2*(2*lp['poly']['dummies'])*lp['poly']['pitch']+66, "H", "M2")
                self.d0.update(101, -65)
                self.d0.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k, -130), 125, 4*lp['poly']['dummies']*lp['poly']['pitch'], self.m1.fingers, "M1")
                add_via_array(cell_i, (101+k, -125), 4*lp['poly']['dummies']*lp['poly']['pitch'], self.m1.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16+k, -125), (self.m1.fingers-1)*2*(2*lp['poly']['dummies'])*lp['poly']['pitch']+66, "H", "M2")
                self.d1.update(101+k, -130)
                self.d1.transform_p(16, 16)

                # 3) gate connection
                add_metal(cell_i, P(0, y_g), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")

                # 4) source to bulk
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                x_len = ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']
                # add_metal(cell_i, P((self.total_fingers-2)*lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), (2*lp['poly']['dummies']+3)*lp['poly']['pitch'], "H", "M2")

                add_metal(cell_i, P(x_len/2, y_s), (self.total_fins+1)*lp["fins"]["pitch"]+205, "V", "M2")

                # 5) d0 to gate connection
                add_v_metal_array(cell_i, P(101, y_g-66), 98, 2*k, self.m0.fingers, "M1")  # body drain contact m0

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))

                self.s.transform_p(x_, y_)
                self.d0.transform_p(x_, y_)
                self.d1.transform_p(x_, y_)

                add_label(self.cell, labels[0], self.d0, "M1LABEL")
                add_label(self.cell, labels[1], self.d1, "M1LABEL")
                add_label(self.cell, labels[2], self.s, "M1LABEL")

            elif pattern == 5: # Common-Centroid -----------------------------------------
                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                poly_no = self.total_fingers*2*(lp['poly']['dummies'])+1
                k = lp['poly']['dummies']*lp['poly']['pitch'] # gap between the same mosfets
                k2 = self.m0.fingers*lp['poly']['dummies']*lp['poly']['pitch']
                k3 = k2 + self.m1.fingers*2*lp['poly']['dummies']*lp['poly']['pitch']
                for i in range(self.total_fingers):
                    mi = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'multiplier': 1, 'mos_type': self.m0.mos_type})
                    mos_i = create_mos(mi, fabric_on = False, con=[0,1,1,1])
                    add_transformed_polygons(mos_i, cell_i,(i*(2*lp['poly']['dummies']) * lp['poly']['pitch'] , 0))

                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions

                # 1) source connection
                add_metal(cell_i, P(0, y_s), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                self.s.update(23 + 16, y_s + 16)

                # 2) drain connection
                add_v_metal_array(cell_i, P(101, -65), 60, 2*k, self.m0.fingers//2, "M1")
                add_via_array(cell_i, (101, -60), 2*k, self.m0.fingers//2, 'H', "V1")
                self.d0.update(101, -65)
                self.d0.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k2, -130), 125, 2*k, self.m1.fingers, "M1")
                add_via_array(cell_i, (101+k2, -125), 2*k, self.m1.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16+k2, -125), ((self.m1.fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']+66, "H", "M2")
                self.d1.update(101+k2, -130)
                self.d1.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k3, -65), 60, 2*k, self.m0.fingers//2, "M1")
                add_via_array(cell_i, (101+k3, -60), 2*k, self.m0.fingers//2, 'H', "V1")

                add_metal(cell_i, P(101-16, -60), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']+66, "H", "M2")

                # 3) gate connection
                add_metal(cell_i, P(0, y_g), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")

                # 4) source to bulk
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                x_len = ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']+66
                add_metal(cell_i, P((self.total_fingers-2)*lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), (2*lp['poly']['dummies']+3)*lp['poly']['pitch'], "H", "M2")

                v_metals = math.ceil(self.total_fingers/8)
                for i in range(v_metals):
                    delta_x = x_len/(v_metals+1)
                    add_metal(cell_i, P((i+1)*delta_x, y_s), (self.total_fins+1)*lp["fins"]["pitch"]+205, "V", "M2")

                # 5) d0 to gate connection
                add_v_metal_array(cell_i, P(101, y_g-66), 98, 2*k, self.m0.fingers//2, "M1")  # body drain contact m0
                add_v_metal_array(cell_i, P(101+k3, y_g-66), 98, 2*k, self.m0.fingers//2, "M1")  # body drain contact m0

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))


                self.s.transform_p(x_, y_)
                self.d0.transform_p(x_, y_)
                self.d1.transform_p(x_, y_)

                add_label(self.cell, labels[0], self.d0, "M1LABEL")
                add_label(self.cell, labels[1], self.d1, "M1LABEL")
                add_label(self.cell, labels[2], self.s, "M1LABEL")