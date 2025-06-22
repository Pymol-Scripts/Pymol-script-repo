'''
See more here: http://www.pymolwiki.org/index.php/stereo_ray
'''

from pymol import cmd


def stereo_ray(filename, width=0, height=0, quiet=1):
    '''
 DESCRIPTION

    "stereo_ray" ray-traces the current scene twice (separated by 
    a six-degree rotation around the y axis)
    and saves a pair of images that can be combined in any image
    manipulation software to form a stereoimage.
    The first argument, the output file name, is mandatory.
    The second and third arguments, the size of the image, are not.
    If the width is given, the height will be calculated.

 USAGE

    stereo_ray filename [, width [, height]]

 EXAMPLE

    stereo_ray output, 1000, 600
    stereo_ray secondImage.png
    '''
    if filename.lower().endswith('.png'):
        filename = filename[:-4]

    cmd.ray(width, height, angle=-3)
    cmd.png(filename + "_r", quiet=quiet)
    cmd.ray(width, height, angle=3)
    cmd.png(filename + "_l", quiet=quiet)

cmd.extend('stereo_ray', stereo_ray)

# vi:expandtab:sw=3
