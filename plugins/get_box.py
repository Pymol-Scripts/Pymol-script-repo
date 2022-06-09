from pymol import cmd
import numpy as np


@cmd.extend
def get_box(box_sel, box_margin):
    """
Create a box (useful for VINA or DOCKTHOR).
USAGE:
    get_box sel, margin
EXAMPLE:
    get_box resn CU and chain A, 5
    
    """
    box_margin = int(box_margin)
    
    box_coords = cmd.get_coords(box_sel)

    max = np.max(box_coords, axis=0) + box_margin
    min = np.min(box_coords, axis=0) - box_margin

    half_size = (max - min) / 2
    center = min + half_size

    size_x, size_y, size_z = half_size * 2
    center_x, center_y, center_z = center

    size_x, size_y, size_z = (
        round(float(size_x), 2),
        round(float(size_y), 2),
        round(float(size_z), 2),
    )

    center_x, center_y, center_z = (
        round(float(center_x), 2),
        round(float(center_y), 2),
        round(float(center_z), 2),
    )
    print("Size:", size_x, size_y, size_z)
    print("Center:", center_x, center_y, center_z)
