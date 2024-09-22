import os

import gdstk
import yaml
import matplotlib.pyplot as plt
from utilities.pvector import P
# from core.constants import Const as C
from shapely.geometry import box
from shapely.geometry import Point, Polygon
import numpy as np

base_dir = os.path.dirname(os.path.abspath(__file__))
f = n_file = os.path.join(base_dir, '../layer_param.yml')

with open(f, 'r') as file:
    lp = yaml.safe_load(file)  # layers properties


def write_gds(cell, name):
    lib = gdstk.Library(unit=1e-9, precision=1e-12)
    lib.add(cell)
    lib.write_gds("outputs/"+name+".gds")


def add_region(cell, bottom_left_coord: tuple, top_right_coord: tuple, name: str):
    x_i, y_i = bottom_left_coord[0], bottom_left_coord[1]
    x_f, y_f = top_right_coord[0], top_right_coord[1]
    layer_no = lp[name]['layer_no'].split(':')
    region = gdstk.rectangle((x_i, y_i), (x_f, y_f), layer=int(layer_no[0]), datatype=int(layer_no[1]))
    cell.add(region)

def add_region_array(cell, bottom_left_coord, dimension, pitch, number, direction, name):
    x_i = bottom_left_coord[0]
    y_i = bottom_left_coord[1]

    for i in range(number):
        add_region(cell, (x_i, y_i),
                   (x_i + dimension[0], y_i + dimension[1]), name)

        if direction == "H":
            x_i += pitch
        elif direction == "V":
            y_i += i*pitch


def add_metal(cell, pos, length, direction, name, return_rect=False):
    m_thickness = lp[name]['thickness']
    rect = 0
    if direction == "H":
        r = [[pos.x, pos.y], [pos.x + length, pos.y+m_thickness]]
    else:
        r = [[pos.x, pos.y], [pos.x + m_thickness, pos.y + length]]

    add_region(cell, (r[0][0], r[0][1]), (r[1][0], r[1][1]), name)

    if return_rect:
        return r


def add_v_metal_array(cell,pos, length: float, pitch:float, number: int, name: str):
    x_i = pos.x
    y_i = pos.y

    for i in range(number):
        add_metal(cell, P(x_i, y_i),length,"V", name)
        x_i += pitch

def go_up(self, cell, pos, bottom_layer, top_layer):
    ...


def go_down(self, cell, pos, top_layer, bottom_layer):
    ...


def add_fabric(cell, poly_no, fin_no):
    poly_height = fin_no * lp['fins']['pitch']
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


def add_power_grid(cell, loc_i, loc_j, m_h, m_v, h_pitch, v_pitch):
    ver_metals = int((loc_j.x - loc_i.x)/h_pitch)
    hor_metals = int((loc_j.y - loc_i.y)/v_pitch)

    for i in range(ver_metals):
        add_metal(cell, P(loc_i.x + i*h_pitch, loc_i.y), hor_metals*v_pitch, "V", m_h)

    for i in range(hor_metals):
        add_metal(cell, P(loc_i.x, loc_i.y + i*h_pitch), ver_metals*v_pitch, "H", m_v)


def add_transformed_polygons(source_cell, target_cell, offset):
    for polygon in source_cell.polygons:
        transformed_points = [(x + offset[0], y + offset[1]) for x, y in polygon.points]
        transformed_polygon = gdstk.Polygon(transformed_points, layer=polygon.layer, datatype=polygon.datatype)
        target_cell.add(transformed_polygon)

    for label in source_cell.labels:
        new_position = (label.origin[0] + offset[0], label.origin[1] + offset[1])
        new_label = gdstk.Label(
            label.text,
            new_position,
            layer=label.layer,
            texttype=label.texttype,
            anchor=label.anchor,
            rotation=label.rotation,
            magnification=label.magnification,
            x_reflection=label.x_reflection
        )
        target_cell.add(new_label)
def add_label(cell, label: str, position: P, name: str):
    layer_no = lp[name]['layer_no'].split(':')
    text_label = gdstk.Label(label, (position.x, position.y), layer=int(layer_no[0]), texttype=int(layer_no[1]))
    cell.add(text_label)


def show_layout(cell, fig_size=(6.4, 4.8)):
    # Visualize the layout using matplotlib
    layer_colors = lp['Colors']
    fig, ax = plt.subplots(figsize=fig_size, dpi=300)

    for polygon in cell.polygons:
        points = polygon.points
        layer = polygon.layer
        color = layer_colors.get(layer, 'gray')  # Default color is black if layer not specified
        ax.fill(*zip(*points), color=color, alpha=0.3)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()

def connect(cell, reg1, reg2, direction, metal_i, metal_j):
    # reg1.x_i  ---  reg2.x_i
    # reg1.y_i  ---  reg2.y_i
    # reg1.x_f  ---  reg2.x_f
    # reg1.x_f  ---  reg2.x_f
    # def add_metal(cell, pos, length, direction, name):

    if(direction == "H"):
        add_metal(cell, P(reg1.x_f, reg1.y_f), reg2.x_i - reg1.x_f, "H", metal_i)

    else:
        print("This has not implemented yet")

def generate_pattern(mos1,mos2,pat):
    pattern = []
    m1 = mos1.multiplier
    m2 = mos1.multiplier
    if pat == "I":
        for i in range(m1):
            pattern.append(mos1)
            pattern.append(mos2)

    elif pat == "N":
        for i in range(m1):
            pattern.append(mos1)
        for i in range(m2):
            pattern.append(mos2)

    return pattern


def pattern_generator(total_fins):
    fingers = []
    for i in range(2, total_fins // 2 + 1):
        if total_fins % i == 0:
            if i % 2 == 0:
                fingers.append(i)

    return fingers


def create_empty_cell(name: str, unit: float, precision: float):
    lib = gdstk.Library(name=name, unit=unit, precision=precision)
    cell = lib.new_cell(name)
    return cell


def fill_v_metals(cell, rect: list, name: str):
    bottom_left_coord, top_right_coord = rect[0], rect[1]

    x_sp = top_right_coord[0] - bottom_left_coord[0] - 2*20
    y_sp = top_right_coord[1] - bottom_left_coord[1]

    via_sp = lp["V0"]["length"] + lp["V0"]["pitch"]

    n_x = max((x_sp//via_sp), 1)
    n_y = max((y_sp//via_sp), 1)

    x_i = bottom_left_coord[0] + 20 + (x_sp - n_x*via_sp + lp["V0"]["pitch"])/2
    y_i = bottom_left_coord[1] + y_sp/2

    for i in range(int(n_x)):
        add_metal(cell, P(x_i + i*via_sp, y_i - 40), 80, "V", name)


def fill_with_v_M1_array(cell, rect: list, pitch, n, name: str):
   for i in range(n):
       fill_v_metals(cell, [(rect[0][0]+i*pitch, rect[0][1]), (rect[1][0]+i*pitch, rect[1][1])], name)

def get_intersection(rect1, rect2):
    # Convert rectangles to Shapely box objects
    rect1_box = box(*rect1[0], *rect1[1])
    rect2_box = box(*rect2[0], *rect2[1])

    # Compute intersection
    intersection = rect1_box.intersection(rect2_box)

    return intersection

class Contact_region:
    def __init__(self, id, rect, metal="x"):
        self.id = id
        self.rect = rect
        self.metal = metal

    def translate(self, xoff, yoff):
        rect_ = [[0,0],[0,0]]
        rect_[0][0] = self.rect[0][0] + xoff
        rect_[1][0] = self.rect[1][0] + xoff
        rect_[0][1] = self.rect[0][1] + yoff
        rect_[1][1] = self.rect[1][1] + yoff
        self.rect = rect_

def find_center_n_span_rects(rectangles):
    ''' Just give a list of rectangles with two corners and get
    the center point both in y and x in the rectangle '''
    def find_max_min(rectangles):
        x_min = float('inf')
        x_max = float('-inf')
        y_min = float('inf')
        y_max = float('-inf')
        for r in rectangles:
            if min(r[0][0], r[1][0]) < x_min:
                x_min = min(r[0][0], r[1][0])
            if max(r[0][0], r[1][0]) > x_max:
                x_max =  max(r[0][0], r[1][0])

            if min(r[0][1], r[1][1]) < y_min:
                y_min = min(r[0][1], r[1][1])
            if max(r[0][1], r[1][1]) > y_max:
                y_max = max(r[0][1], r[1][1])

        return [[x_min, y_min], [x_max, y_max]]

    def rectangle_center(rect):
        center_x = (rect[1][0] + rect[0][0]) / 2
        center_y = (rect[1][1] + rect[0][1]) / 2
        return (center_x, center_y)

    def is_point_inside_rect(rect, pt):
        rectangle = Polygon([rect[0], (rect[1][0], rect[0][1]), rect[1], (rect[0][0], rect[1][1])])
        point = Point(pt[0], pt[1])
        return rectangle.contains(point)

    best_center = 0
    centers = [rectangle_center(rect) for rect in rectangles]
    median_x = np.median([center[0] for center in centers])
    median_y = np.median([center[1] for center in centers])
    median_center = (median_x, median_y)

    for rect in rectangles:
        if is_point_inside_rect(rect, median_center):
            best_center = median_center
            break

        else:
            best_center = centers[np.argmin(np.linalg.norm(np.array(centers)-np.array(median_center),  axis=1))]

    return best_center, find_max_min(rectangles)



def divide_line(xi:int, xf:int, n: int, ext: float, m: str):
    '''
    xi: initial point, xf: final point n: number of metals ext: extension m: metal string
    return: list, the initial coordinates on the line
    '''
    l = xf - xi
    width = lp[m]["width"]
    pitch = lp[m]["pitch"]
    delta = l - 2*ext
    num = min((delta-width)//pitch + 1, n)

    xi = (delta - ((num-1)*pitch + width))/2 + ext + xi
    li = [xi+i*pitch for i in range(num)]

    return li


def m1_array(xi, xf, stack=1, n = 1):
    delta = xf-xi
    n = ((delta - lp["M1"]["width"])//(stack*lp["M1"]["pitch"])) + 1
    arr = [xi+i*lp["M1"]["pitch"] for i in range(n)]
    return arr


class Mos:
    def __init__(self, p: dict):
        self.fins = p["fins"]
        self.fingers = p["fingers"]
        self.stack = p["stack"]
        self.multiplier = p["multiplier"]
        self.id = p["id"]
        self.mos_type = p["mos_type"]
