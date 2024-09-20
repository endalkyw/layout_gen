import numpy as np
import pandas as pd
import shapely
from utilities.elements_utils import divide_line
import matplotlib.pyplot as plt

def main():
    p0 = [-202.0, 3634.0]
    p1 = [40, 40]

    l = divide_line(p0[0], p0[1], 20, 20, "M2")
    print(l)

    plt.plot(p0, p1)
    for li in l:
        plt.plot([l],[40], 'ro')

    plt.show()

    print(l)


if __name__ == "__main__":
    main()