'''
http://pymolwiki.org/index.php/load_img_stack

PyMOL image stack loading module

(c) Thomas Holder, Schrodinger, Inc.
'''

from pymol import cmd, CmdException

def imgframeiter(img):
    '''
    Iterate over frames of a multi-frame image (like tiff)
    '''
    import itertools
    try:
        for i in itertools.count():
            img.seek(i)
            yield img
    except EOFError:
        raise StopIteration

def load_img_stack(pattern, name='', grid=1.0, channel=0, normalize=1,
        extent=None, quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a stack of images as a map

ARGUMENTS

    pattern = str: image filename or pattern

    name = str: map object name to create

    grid = float: grid spacing in Angstrom {default: 1.0}

    channel = int: color channel for RGB images {default: 0}

    normalize = 0 or 1: normalize data {default: 1}

    extent = 3-float: (a,b,c) edge lengths in Angstrom, overwrites "grid"
    arguments if given {default: }

EXAMPLES

    load_img_stack img*.tif, extent=(21.0, 14.5, 18.2)
    '''
    import glob
    import numpy
    from chempy.brick import Brick

    try:
        from PIL import Image
    except ImportError:
        import Image

    channel, normalize, quiet = int(channel), int(normalize), int(quiet)

    if not name:
        name = _self.get_unused_name('map')

    if isinstance(grid, str):
        grid = _self.safe_eval(grid)
    if not isinstance(grid, (tuple, list)):
        grid = (grid,) * 3

    stack = []
    size = None

    filenames = glob.glob(_self.exp_path(pattern))

    if not filenames:
        raise CmdException('no such files')

    for filename in sorted(filenames):
        img = Image.open(filename)
        if size is None:
            size = img.size
        for img in imgframeiter(img):
            if img.size != size:
                if not quiet:
                    print('Image size mismatch: %s != %s' % (img.size, size))
                continue
            a = numpy.reshape(img, (size[0], size[1], -1))
            stack.append(a[..., channel])

    stack = numpy.asfarray(stack)
    stack = stack.swapaxes(0, 2)[:, ::-1, ::-1]

    if min(stack.shape) < 2:
        raise CmdException('insufficient grid dimensions: ' + str(stack.shape))

    if normalize:
        stack -= stack.mean()
        stack /= stack.std()

    if extent:
        if isinstance(extent, str):
            extent = _self.safe_eval(extent)
        grid = [(float(e) / (s - 1)) for (e, s) in zip(extent, stack.shape)]
        if not quiet:
            print(' Setting grid = ' + str(grid))

    brick = Brick.from_numpy(stack, grid)
    _self.load_brick(brick, name)

cmd.extend('load_img_stack', load_img_stack)
