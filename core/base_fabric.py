from utilities.elements_utils import *
from utilities.pvector import P
def create_fabric(cell, poly_no, fin_no):
    poly_height = fin_no * lp['fins']['pitch']-lp['fins']['pitch']+lp['fins']['width']
    fin_len = (poly_no - 1) * lp['poly']['pitch']

    for i in range(fin_no):
        x_initial = 23
        y_init = i * lp['fins']['pitch']
        add_region(cell, (x_initial, y_init), (x_initial + fin_len, y_init + lp['fins']['width']), "RXFIN")

    for i in range(poly_no):
        y_init = 0
        x_init = x_initial + i * lp['poly']['pitch']
        add_region(cell, (x_init - 23, y_init), (x_init + lp['poly']['width'] + 23, y_init + poly_height), "TI")
        add_region(cell, (x_init, y_init), (x_init + lp['poly']['width'], y_init + poly_height), "poly")


def get_coordinate(poly_number, fin_number):
    return P(poly_number*lp["poly"]["pitch"] + 23, fin_number * lp["fins"]["pitch"])




class ConnectionRects:
    def __init__(self):
        self.G = []
        self.B = []
        self.D = []
        self.S = []
        self.dim = []


