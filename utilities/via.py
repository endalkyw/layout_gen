import gdstk
import yaml
from utilities.pvector import P
from core.constants import Const as C
from utilities.elements_utils import *

f = 'layer_param.yml'
with open(f, 'r') as file:
    lp = yaml.safe_load(file)  # layers properties


def add_via(cell, bottom_left_coord: tuple, name: str):
    x_i, y_i = bottom_left_coord[0], bottom_left_coord[1]
    l = lp[name]["length"]
    layer_no = lp[name]['layer_no'].split(':')
    via_i = gdstk.rectangle((x_i, y_i), (x_i + l, y_i + l), layer=int(layer_no[0]), datatype=int(layer_no[1]))
    cell.add(via_i)


def add_via_array(cell, bottom_left_crd: tuple, pitch: float, number: int, dir: str, name: str):
    x_i = bottom_left_crd[0]
    y_i = bottom_left_crd[1]

    for i in range(number):
        add_via(cell, (x_i, y_i,), name)
        if dir == 'H':
            x_i += pitch
        else:
            y_i += pitch


def fill_with_vias(cell, di: str, bottom_left_coord: tuple, top_right_coord: tuple, name: str):
    h_margin = lp[name]["length"]
    v_margin = lp[name]["length"]

    if di == "H":
        available_h_space = top_right_coord[0]-bottom_left_coord[0] - 2*h_margin
        num_vias = int(available_h_space/126) + 1
        x_gap = available_h_space - (num_vias-1)*126 - lp[name]['length']
        available_v_space = top_right_coord[1] - bottom_left_coord[1]
        y_gap = available_v_space - lp[name]['length']
        for i in range(num_vias):
            add_via(cell, (bottom_left_coord[0] + i*126 + x_gap/2 + h_margin, bottom_left_coord[1]+y_gap/2), name)
    else:
        available_v_space = top_right_coord[1] - bottom_left_coord[1] - 2 * v_margin
        num_vias = int(available_v_space / 126) + 1
        y_gap = available_v_space - (num_vias - 1) * 126 - lp[name]['length']
        available_h_space = top_right_coord[0] - bottom_left_coord[0]
        x_gap = available_h_space - lp[name]['length']
        for i in range(num_vias):
            add_via(cell, (bottom_left_coord[0] + x_gap / 2, bottom_left_coord[1] + y_gap / 2 + i * 126 + x_gap / 2 + h_margin), name)


def fill_area_vias(cell, rect: list, name: str):
    bottom_left_coord, top_right_coord = rect[0], rect[1]

    x_sp = top_right_coord[0] - bottom_left_coord[0] - 2*C.v_mrgn
    y_sp = top_right_coord[1] - bottom_left_coord[1] - 2*C.v_mrgn
    via_sp = lp[name]["length"] + lp[name]["pitch"]

    n_x = max((x_sp//via_sp), 1)
    n_y = max((y_sp//via_sp), 1)
    x_i = bottom_left_coord[0] + C.v_mrgn + (x_sp - n_x*via_sp + lp[name]["pitch"])/2
    y_i = bottom_left_coord[1] + C.v_mrgn + (y_sp - n_y*via_sp + lp[name]["pitch"])/2

    for i in range(int(n_x)):
        for j in range(int(n_y)):
            add_via(cell, (x_i + i*via_sp, y_i + j*via_sp), name)


def add_vias(cell, di: str, bottom_left_coord: tuple, top_right_coord: tuple, name: str):
    h_margin = lp[name]["length"]
    v_margin = lp[name]["length"]

    if di == "H":
        available_h_space = top_right_coord[0]-bottom_left_coord[0] - 2*h_margin
        num_vias = int(available_h_space/126) + 1
        x_gap = available_h_space - (num_vias-1)*126 - lp[name]['length']
        available_v_space = top_right_coord[1] - bottom_left_coord[1]
        y_gap = available_v_space - lp[name]['length']
        for i in range(num_vias):
            add_via(cell, (bottom_left_coord[0] + i*126 + x_gap/2 + h_margin, bottom_left_coord[1]+y_gap/2), name)
    else:
        available_v_space = top_right_coord[1] - bottom_left_coord[1] - 2 * v_margin
        num_vias = int(available_v_space / 126) + 1
        y_gap = available_v_space - (num_vias - 1) * 126 - lp[name]['length']
        available_h_space = top_right_coord[0] - bottom_left_coord[0]
        x_gap = available_h_space - lp[name]['length']
        for i in range(num_vias):
            add_via(cell, (bottom_left_coord[0] + x_gap / 2, bottom_left_coord[1] + y_gap / 2 + i * 126 + x_gap / 2 + h_margin), name)


def fill_with_vias_array(cell, di: str, bottom_left_coord: tuple, top_right_coord: tuple, n: int, pitch: float, name: str):
    for i in range(n):
        r = [(bottom_left_coord[0]+i*pitch, bottom_left_coord[1]), (top_right_coord[0]+i*pitch, top_right_coord[1])]
        fill_area_vias(cell, r, name)


def add_vias_at_intersection(cell, rect, metal):
    l = lp[metal]["layer_no"].split(":")
    polygon_1 = gdstk.Polygon([rect[0], (rect[1][0], rect[0][1]), rect[1],(rect[0][0],rect[1][1])])

    for polygon_2 in cell.polygons:
        if polygon_2.layer == int(l[0]) and polygon_2.datatype == int(l[1]):
            intersection_polygons = gdstk.boolean([polygon_1], [polygon_2], "and")
            if intersection_polygons:
                p = intersection_polygons[0].points
                if metal == "M1":
                    add_via(cell, p[2],"V1")
                if metal == "M2":
                    add_via(cell, p[2],"V2")

                # for debugging purpose
                # rect = gdstk.rectangle((p[0][0], p[0][1]),(p[2][0],p[2][1]))
                # cell.add(rect)


def add_vias_at_recti_v_rectj(cell, rect_i, rect_j, metal):
    polygon_1 = gdstk.Polygon([rect_i[0], (rect_i[1][0], rect_i[0][1]), rect_i[1],(rect_i[0][0],rect_i[1][1])])
    polygon_2 = gdstk.Polygon([rect_j[0], (rect_j[1][0], rect_j[0][1]), rect_j[1],(rect_j[0][0],rect_j[1][1])])
    intersection_polygons = gdstk.boolean([polygon_1], [polygon_2], "and")

    if intersection_polygons:
        p = intersection_polygons[0].points
        if metal == "M1":
            add_via(cell, p[2],"V1")
        elif metal == "M2":
            add_via(cell, p[2],"V2")
        elif metal == "M3":
            add_via(cell, p[2],"V3")