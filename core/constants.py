# constants.py

class Const:
    # general ----------------------------
    m_v_ext = 20 # metal extension from the via edge
    m_m_ext = 20 # H and V metal connection extension
    v_mrgn  = 20 # via margin, metal edge to via distance
    mm_overlap = 10 # overlap length between two same layer metals

    # for the base_cell generator
    dist_0  = 44
    gate_h0 = 80
    gate_h1 = 160
    rx_m1   = 23
    rx_po   = 71
    bulk_h1 = 284  # bulk height yf + bulk_h
    bulk_h2 = 380  #320 the bulk height yf + bulk_h
    po_cb   = 20   # cb extension from the cb
    ct_h0  = 162  # 115 107 poly cut height from y_f
    ct_h1  = 258  # 258 200, 180 poly cut height from y_f 2
    ct_hb  = 80   # bottom poly cut dist from y_i

    gate_m1_len = 80  # m1 connector gate length

    v_mult_space   = 7
    v_mult_space_2 = 9
    h_mult_space   = 3
    # mos.py - to generate a single mosfet

