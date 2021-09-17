import numpy as np
import sys

if __name__ == "__main__":
    n_points = int(sys.argv[1])
    points = np.random.uniform(0, 1, size = (n_points, 3))

    with open("points.txt", "w") as f:
        [f.write(f"{point}\n") for point in points]