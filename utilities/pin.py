import numpy as np
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt

def rectangle_center(rect):
    center_x = (rect[1][0]+rect[0][0])/2
    center_y = (rect[1][1]+rect[0][1])/2
    return (center_x, center_y)


def is_point_inside_rect(rect, pt):
    rectangle = Polygon([rect[0], (rect[1][0], rect[0][1]), rect[1], (rect[0][0], rect[1][1])])
    point = Point(pt[0], pt[1])
    return rectangle.contains(point)

def find_best_point(rectangles):
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

    return best_center, median_center, centers

def show_rectangles(rectangles, best_center, median_center, centers):
    # Step 4: Plotting the rectangles and the best center
    fig, ax = plt.subplots()

    # Plot each rectangle
    for corner1, corner2 in rectangles:
        rect_polygon = Polygon([
            (corner1[0], corner1[1]),  # bottom-left
            (corner2[0], corner1[1]),  # bottom-right
            (corner2[0], corner2[1]),  # top-right
            (corner1[0], corner2[1])  # top-left
        ])
        x, y = rect_polygon.exterior.xy
        ax.fill(x, y, alpha=0.4, label=f"Rectangle ({corner1}, {corner2})")

    # Plot the best center
    ax.plot(best_center[0], best_center[1], 'ro', label="Best Center", markersize=10)
    ax.plot(median_center[0], median_center[1], 'r*', label="Best Center", markersize=10)

    for c in centers:
        ax.plot(c[0], c[1], 'bo', label="Best Center", markersize=5)

    # Add labels and legend
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    # ax.legend(loc="best")

    # Show the plot
    plt.title("Rectangles and Best Center")
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()


def main():
    rectangles = []
    for i in range(5):
        rectangles.append(np.random.randint(100)*np.random.rand(2, 2))

    best_centers, median_center, centers = find_best_point(rectangles)
    show_rectangles(rectangles, best_centers,median_center, centers)


if __name__ == "__main__":
    main()