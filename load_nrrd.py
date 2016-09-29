'''
PyMOL NRRD map file loading

(c) Thomas Holder, Schrodinger, Inc.

License: BSD-2
'''

from pymol import cmd, CmdException

@cmd.extend
def load_nrrd(filename, prefix='channel', channel='', _self=cmd):
    '''
DESCRIPTION

    Load maps in NRRD file format.

    Spec:
    http://teem.sourceforge.net/nrrd/format.html

    Example files:
    http://www.cs.utah.edu/~jmk/simian/download.htm
    '''
    import numpy
    import shlex
    from chempy.brick import Brick

    handle = open(filename, 'rb')

    magic = handle.readline()
    assert magic.startswith('NRRD')

    header = {}

    while True:
        line = handle.readline()
        if not line.strip():
            break
        if line.startswith('#'):
            continue
        key, value = line.split(':')
        header[key] = value.strip()

    assert header['encoding'] == 'raw'

    endian = '>' if header.get('endian') == 'big' else '<'
    dtype = {
        'signed char': 'i1',
        'int8': 'i1',
        'int8_t': 'i1',
        'uchar': 'u1',
        'unsigned char': 'u1',
        'uint8': 'u1',
        'uint8_t': 'u1',
        'float': 'f4',
        'double': 'f8',
    }[header['type']]

    data = numpy.fromfile(handle, endian + dtype)
    handle.close()

    sizes = [int(i) for i in header['sizes'].split()]
    grid = [float(i) for i in header.get('spacings', '0 1 1 1').split()[1:]]
    channels = shlex.split(header.get('labels', 'abcdefg'))[0]

    assert len(sizes) == 4

    sizes = tuple(reversed(sizes))
    grid = tuple(reversed(grid))

    data = data.reshape(sizes)

    for i in range(sizes[-1]):
        ch = channels[i]
        name = prefix

        if not channel:
            name += '_' + ch
        elif channel != ch:
            continue

        stack = data[..., i].astype('f4')

        # normalize
        stack -= stack.mean()
        stack /= stack.std()

        brick = Brick.from_numpy(stack, grid)
        _self.load_brick(brick, name)

if __name__ == '__main__':
    # pymol -r load_nrrd.py -- file.nrrd
    import sys
    load_nrrd(sys.argv[1])
