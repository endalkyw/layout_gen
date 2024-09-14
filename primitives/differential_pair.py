from utilities.elements_utils import *
from core import create_base_layout, create_fabric
from utilities.pvector import P
from utilities.via import *
from primitives.mos import *

from core.constants import Const as C
class differential_pair:
    def __init__(self, m0: Mos, m1: Mos, name):
        assert m0.fins >= 2, "Number of fins must be greater than 1"
        assert m0.fingers == m1.fingers, "Fingers must be the same"
        assert m0.fins == m1.fins, "fins must be the same"
        assert m0.mos_type == m1.mos_type, "MOS type must be the same"

        self.m0 = m0
        self.m1 = m1
        self.total_fingers = m0.fingers + m1.fingers
        self.stack = m0.stack
        self.total_fins = m0.fins
        self.mos_type = m0.mos_type
        self.cell = create_empty_cell(name,  1e-9, 1e-12)

    def pin_polygon(self, rect: list, label: str, metal: str):
        pass

    def add_labels(self):
        pass




    def create_layout(self, pattern: int, labels =["d0", "d1", "g0", "g1", "s", "b"], label_pos =["c", "c", "c", "c", "c"],
                        con =[1, 1, 1, 1, 1], row_number = 1, fabric_on=True):

        fin_no = self.total_fins + 4 * lp['fins']['dummies']
        poly_no = self.stack*self.total_fingers + 2 * lp['poly']['dummies']

        cell_i = create_empty_cell("dp", 1e-9, 1e-12)
        n0 = self.m0.fingers//2 + int(self.m0.fingers % 2)
        k1 = C.rx_m1 + self.stack * lp["M1"]["pitch"]


        if pattern == 0 or pattern == 2 or pattern == 4:  # no diffusion break


            if pattern == 0:
                create_base_layout(cell_i, self.total_fingers, self.total_fins, self.mos_type, self.stack,"dp_0")
                n0 = self.m0.fingers // 2 + int(self.m0.fingers % 2)
                n1 = self.m1.fingers // 2 + int(self.m1.fingers % 2)

                # drain 0 connection
                add_v_metal_array(cell_i, P(k1, -65), 66, self.stack*2 * lp['M2']['pitch'], n0, "M1")
                add_via_array(cell_i, (k1, -60), self.stack*2 * lp['M2']['pitch'], n0, 'H', "V1")
                add_metal(cell_i, P(k1 - 10, -60), self.stack*(n0 - 1) * 2 * lp['M2']['pitch'] + lp['M1']['width'] + 20, "H", "M2")
                d0 = P(101 + 16, -60 + 16)

                # drain 1 connection
                k = self.stack*n0 * 2 * lp['M2']['pitch']
                add_v_metal_array(cell_i, P(k1 + k, -65), 66, self.stack*2 * lp['M1']['pitch'], n1, "M1")
                add_via_array(cell_i, (k1 + k, -60), self.stack*2 * lp['M2']['pitch'], n1, 'H', "V1")
                add_metal(cell_i, P(k1 + k - 10, -60), (n1 - 1) * self.stack*2 * lp['M2']['pitch'] + lp['M1']['width'] + 20, "H", "M2")

                d1 = P(k1 + k + 16, -60 + 16)
                y_g = self.total_fins * lp['fins']['pitch'] + 65
                xi = 62
                len = self.stack*(self.total_fingers / 2) * lp['M2']['pitch'] + 16 - xi

                # gate0 connections
                add_metal(cell_i, P(xi, y_g), len, "H", "M1")
                g0 = P(xi + 16, y_g + 16)

                # gate1 connections
                add_metal(cell_i, P(len + xi + 46, y_g), len, "H", "M1")
                g1 = P(len + xi + 46 + 16, y_g + 16)

                b = P(0 + self.stack*self.total_fingers*lp['poly']['pitch']/2, self.total_fins*lp['fins']['pitch'] + 264)

            elif pattern == 2:
                create_base_layout(cell_i, self.total_fingers, self.total_fins, self.mos_type, self.stack,"dp_1")
                d0, d1, g0, g1 = self.interdigitated(cell_i)
                b = P(0 + self.total_fingers*lp['poly']['pitch']/2, self.total_fins*lp['fins']['pitch'] + 310)

            elif pattern == 4:
                create_base_layout(cell_i, self.total_fingers, self.total_fins, self.mos_type, self.stack,"dp_2")
                d0, d1, g0, g1 = self.commoncentroid(cell_i)
                b = P(0 + self.total_fingers*lp['poly']['pitch']/2, self.total_fins*lp['fins']['pitch'] + 310)

            y_s = 50


            # common to all connections: source, bulk (already connected)
            add_metal(cell_i, P(0, y_s), self.stack*(self.total_fingers+1)*lp['poly']['pitch'], "H", "M2")
            add_via_array(cell_i, (23, y_s), self.stack*2*lp['M1']['pitch'], self.total_fingers//2+1, 'H', "V1")
            s = P(0 + 39, y_s+16)

            n0 = self.m0.fingers//2 + int(self.m0.fingers % 2)
            n1 = self.m1.fingers//2 + int(self.m1.fingers % 2)


            x_ = 0
            y_ = 0
            if fabric_on:
                create_fabric(self.cell, poly_no, fin_no)
                x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width']) / 2

            add_label(cell_i, labels[0], d0, "M1LABEL")
            add_label(cell_i, labels[1], d1, "M1LABEL")
            add_label(cell_i, labels[2], g0, "M1LABEL")
            add_label(cell_i, labels[3], g1, "M1LABEL")
            add_label(cell_i, labels[4], s,  "M1LABEL")
            add_label(cell_i, labels[5], b,  "M1LABEL")
            add_transformed_polygons(cell_i, self.cell, (x_, y_))


        elif pattern == 1 or pattern == 3 or pattern == 5: # diffusion break

            if pattern == 1: # Clustered -------------------------------------------------
                mos_a_cell = create_mos(self.m0, fabric_on = False, labels =  {"d":"d0","g":"g0", "b":"b", "s":0})
                mos_b_cell = create_mos(self.m1, fabric_on = False, labels =  {"d":"d1","g":"g1", "b":0, "s":0})
                poly_no = self.total_fingers + 4 * lp['poly']['dummies']-1

                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                add_transformed_polygons(mos_a_cell, cell_i, (0, 0))
                add_transformed_polygons(mos_b_cell, cell_i, ((2*lp['poly']['dummies'] + self.m0.fingers - 1) * lp['poly']['pitch'] , 0))
                
                s = P(0, 0)
                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions 
                # 1) source connection
                add_metal(cell_i, P(self.m0.fingers*lp['poly']['pitch'], y_s), (2*lp['poly']['dummies'])*lp['poly']['pitch'], "H", "M2")
                s.update(23, y_s)
                s.transform_p(16, 16)
                                
                # 4) source to bulk  
                b = P(0,0)
                x_len = (self.m0.fingers + lp['poly']['dummies'])*lp["poly"]["pitch"]
                b.update(((self.m0.fingers+1)/2)*lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200)
                add_metal(cell_i, P(((self.m0.fingers+1)/2)*lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+200), (((self.total_fingers+2)/2)+4)*lp['poly']['pitch'], "H", "M1")
                b.transform_p(16,16)

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))
                
                s.transform_p(x_, y_)
                b.transform_p(x_, y_)

                add_label(self.cell, labels[4], s, "M1LABEL")
                # add_label(self.cell, labels[5], b, "M1LABEL") this should be the thing I need to use (on m2)
            
            elif pattern == 3: # Interdigitated ------------------------------------------
                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                poly_no = self.total_fingers*2*(lp['poly']['dummies'])+1
                g_type = ["dp_1x0","dp_1x1"]
                for i in range(self.total_fingers):
                    mi = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'multiplier': 1, 'mos_type': self.m0.mos_type})
                    mos_i = create_mos(mi, fabric_on = False, d_con = False, gate_connection_type = g_type[i%2])
                    add_transformed_polygons(mos_i, cell_i,(i*(2*lp['poly']['dummies']) * lp['poly']['pitch'] , 0))

                s  = P(0, 0)
                d0 = P(0, 0)
                d1 = P(0, 0)
                g0 = P(0, 0)
                g1 = P(0, 0)
                b  = P(0, 0) 

                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions 
                
                # 1) source connection
                add_metal(cell_i, P(0, y_s), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                s.update(23, y_s)   
                s.transform_p(16, 16)

                # 2) drain connection 
                k = 2*lp['poly']['dummies']*lp['poly']['pitch'] # gap between the mosfets
                add_v_metal_array(cell_i, P(101, -65), 66, 2*k, self.m0.fingers, "M1")
                add_via_array(cell_i, (101, -60), 2*k, self.m0.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16, -60), (self.m0.fingers-1)*2*(2*lp['poly']['dummies'])*lp['poly']['pitch']+66, "H", "M2")  
                d0.update(101, -65)   
                d0.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k, -130), 125, 4*lp['poly']['dummies']*lp['poly']['pitch'], self.m1.fingers, "M1")
                add_via_array(cell_i, (101+k, -125), 4*lp['poly']['dummies']*lp['poly']['pitch'], self.m1.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16+k, -125), (self.m1.fingers-1)*2*(2*lp['poly']['dummies'])*lp['poly']['pitch']+66, "H", "M2")
                d1.update(101+k, -130)   
                d1.transform_p(16, 16)
                
                # 3) gate
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+17), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+97), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                g0.transform_p(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+17)
                g0.transform_p(16, 16)

                g1.transform_p(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+97)
                g1.transform_p(16, 16)

               
                # 4) bulk  
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+249), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                x_len = ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']
                b.update(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+249)
                b.transform_p(16,16)

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))
                
              
                s.transform_p(x_, y_)
                d0.transform_p(x_, y_)
                d1.transform_p(x_, y_)
                g0.transform_p(x_, y_)
                g1.transform_p(x_, y_)
                b.transform_p(x_, y_)

                add_label(self.cell, labels[0], d0, "M1LABEL")
                add_label(self.cell, labels[1], d1, "M1LABEL")
                add_label(self.cell, labels[2], g0, "M1LABEL")
                add_label(self.cell, labels[3], g1, "M1LABEL")
                add_label(self.cell, labels[4], s, "M1LABEL")
                add_label(self.cell, labels[5], b, "M1LABEL")

            elif pattern == 5: # Common-Centroid -----------------------------------------
                cell_i= create_empty_cell("cm", 1e-9, 1e-12)
                poly_no = self.total_fingers*2*(lp['poly']['dummies'])+1
                k = lp['poly']['dummies']*lp['poly']['pitch'] # gap between the same mosfets
                k2 = self.m0.fingers*lp['poly']['dummies']*lp['poly']['pitch']
                k3 = k2 + self.m1.fingers*2*lp['poly']['dummies']*lp['poly']['pitch']
                
                g_type = ["dp_1x0","dp_1x1"]

                for i in range(self.total_fingers):
                    mi = Mos({'id': 'A', 'fins': self.m0.fins, 'fingers': 1, 'multiplier': 1, 'mos_type': self.m0.mos_type})
                    
                    if i < self.m0.fingers//2 or i>= self.m0.fingers//2 + self.m1.fingers:
                        mos_i = create_mos(mi, fabric_on = False, d_con=False, gate_connection_type=g_type[1])
                    else:
                        mos_i = create_mos(mi, fabric_on = False, d_con=False, gate_connection_type=g_type[0])
                    
                    add_transformed_polygons(mos_i, cell_i,(i*(2*lp['poly']['dummies']) * lp['poly']['pitch'] , 0))
                s  = P(0, 0)
                d0 = P(0, 0)
                d1 = P(0, 0)
                g0 = P(0, 0)
                g1 = P(0, 0)
                b  = P(0, 0) 

                y_s = 50
                y_g = self.total_fins*lp['fins']['pitch'] + 65
                # add all the metal connenctions 
                
                # 1) source connection
                add_metal(cell_i, P(0, y_s), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M2")
                s.update(23 + 16, y_s + 16)   
                
                # 2) drain connection 
                add_v_metal_array(cell_i, P(101, -65), 60, 2*k, self.m0.fingers//2, "M1")
                add_via_array(cell_i, (101, -60), 2*k, self.m0.fingers//2, 'H', "V1")
                d0.update(101, -65)
                d0.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k2, -130), 125, 2*k, self.m1.fingers, "M1")
                add_via_array(cell_i, (101+k2, -125), 2*k, self.m1.fingers, 'H', "V1")
                add_metal(cell_i, P(101-16+k2, -125), ((self.m1.fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']+66, "H", "M2")
                d1.update(101+k2, -130)
                d1.transform_p(16, 16)

                add_v_metal_array(cell_i, P(101+k3, -65), 60, 2*k, self.m0.fingers//2, "M1")
                add_via_array(cell_i, (101+k3, -60), 2*k, self.m0.fingers//2, 'H', "V1")
                
                add_metal(cell_i, P(101-16, -60), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']+66, "H", "M2")

                
                # 3) gate
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+17), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+97), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                g0.transform_p(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+17)
                g0.transform_p(16, 16)

                g1.transform_p(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+97)
                g1.transform_p(16, 16)

               
                # 4) bulk  
                add_metal(cell_i, P(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+249), ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch'], "H", "M1")
                x_len = ((self.total_fingers-1)*(2*(lp['poly']['dummies'])))*lp['poly']['pitch']
                b.update(lp['poly']['pitch'], (self.total_fins+1)*lp["fins"]["pitch"]+249)
                b.transform_p(16,16)

                x_ = 0
                y_ = 0
                if fabric_on:
                    create_fabric(self.cell, poly_no, fin_no)
                    x_ = (0 + 1) * (lp['poly']['dummies'] - 1) * lp['poly']['pitch'] + 23 + lp['poly']['width'] / 2
                    y_ = (0 + 1) * lp['fins']['dummies'] * lp["fins"]["pitch"] - (lp['fins']['pitch'] - lp['fins']['width'])/2

                add_transformed_polygons(cell_i, self.cell, (x_, y_))
                
              
                s.transform_p(x_, y_)
                d0.transform_p(x_, y_)
                d1.transform_p(x_, y_)
                g0.transform_p(x_, y_)
                g1.transform_p(x_, y_)
                b.transform_p(x_, y_)

                add_label(self.cell, labels[0], d0, "M1LABEL")
                add_label(self.cell, labels[1], d1, "M1LABEL")
                add_label(self.cell, labels[2], g0, "M1LABEL")
                add_label(self.cell, labels[3], g1, "M1LABEL")
                add_label(self.cell, labels[4], s, "M1LABEL")
                add_label(self.cell, labels[5], b, "M1LABEL")
               





    def interdigitated(self, cell_i):
        n0 = self.m0.fingers // 2 + int(self.m0.fingers % 2)
        n1 = self.m1.fingers // 2 + int(self.m1.fingers % 2)

        nx = self.m0.fingers//2
        # gate0 connection
        xi = 62
        y_g = self.total_fins * lp['fins']['pitch'] + 65
        len = (self.total_fingers) * lp['M2']['pitch'] + 16 - xi
        add_metal(cell_i, P(xi, y_g), len, "H", "M1")
        g0 = P(xi + 16, y_g + 16)

        y_g = self.total_fins * lp['fins']['pitch'] + 145
        add_metal(cell_i, P(xi, y_g), len, "H", "M1")
        g1 = P(xi + 16, y_g + 16)

        k = lp['M2']['pitch']
        # drain 0 connection
        add_v_metal_array(cell_i, P(101, -65), 60, 4*lp['M2']['pitch'], n0, "M1")
        add_via_array(cell_i, (101, -60), 4*lp['M2']['pitch'], n0, 'H', "V1")
        add_metal(cell_i, P(101 - 10, -60), (n0-1)*4*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
        d0 = P(101+16, -60+16)

        # drain 1 connection
        add_v_metal_array(cell_i, P(101+2*k, -140), 135, 4*lp['M1']['pitch'], n1, "M1")
        add_via_array(cell_i, (101+2*k, -135), 4*lp['M2']['pitch'], n1, 'H', "V1")
        add_metal(cell_i, P(101+2*k - 10, -135), (n1-1)*4*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
        d1 = P(101+2*k+16, -135+16)

        return d0, d1, g0, g1

    def commoncentroid(self, cell_i):
        n0 = self.m0.fingers // 2 + int(self.m0.fingers % 2)
        n1 = self.m1.fingers // 2 + int(self.m1.fingers % 2)
        d0 = P(0,0)
        d1 = P(0,0)
        g0 = P(0,0)
        g1 = P(0,0)

        y_g = self.total_fins*lp['fins']['pitch'] + 65
        
        # drain 0_1 connection
        add_v_metal_array(cell_i, P(101, -65), 60, 2*lp['M2']['pitch'], n0//2, "M1")
        add_via_array(cell_i, (101, -60), 2*lp['M2']['pitch'], n0//2, 'H', "V1")
        d0.update(101+16, -60+16)
                
        # drain 1__ connection
        k1 = (n0/2)*2*lp['M2']['pitch']
        add_v_metal_array(cell_i, P(101+k1, -143), 138, 2*lp['M1']['pitch'], n1, "M1")
        add_via_array(cell_i, (101+k1, -138), 2*lp['M2']['pitch'], n1, 'H', "V1")
        add_metal(cell_i, P(101+k1 - 10, -138), (n1-1)*2*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")
        d1.update(101+k1 + 16, -138 + 16)
                  
        # drain 0_2 connection
        k = n0*3*lp['M2']['pitch']
        add_v_metal_array(cell_i, P(101+k, -65), 60, 2*lp['M2']['pitch'], n0//2, "M1")
        add_via_array(cell_i, (101+k, -60), 2*lp['M2']['pitch'], int(n0/2), 'H', "V1")
        
        # gate drain contact m0
        add_metal(cell_i, P(101 - 10, -60), 2*(n1+n0-1)*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M2")

        # gate 0 contact
        add_metal(cell_i, P(101 - 10, y_g), 2*(n1+n0-1)*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M1")
        g0.update(101 - 10, y_g)
        g0.transform_p(16,16)
        
        # gate 1 contact
        add_metal(cell_i, P(101+k1 - 10, y_g+80), (n1-1)*2*lp['M2']['pitch']+lp['M1']['width']+20, "H", "M1")
        g1.update(101+k1 - 10, y_g+80)
        g1.transform_p(16,16)

     
        return d0, d1, g0, g1