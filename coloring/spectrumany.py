'''
(c) 2010 Thomas Holder
'''
 
from pymol import cmd, stored
 
def spectrumany(expression, color_list, selection='(all)', minimum=None, maximum=None):
    '''
DESCRIPTION
 
    Define a color spectrum with as many color-stops as you like (at least 2).
 
USAGE
 
    spectrumany expression, color_list [, selection [, minimum [, maximum ]]]
 
ARGUMENTS
 
    expression = count, resi, b, q, or pc: respectively, atom count, residue
    index, temperature factor, occupancy, or partial charge {default: count}
 
    color_list = string: Space separated list of colors
 
    ... all other arguments like with `spectrum` command
 
EXAMPLE
 
    spectrumany count, forest green yellow white
    spectrumany b, red yellow white, (polymer), maximum=100.0
 
SEE ALSO
 
    spectrum
    '''
    colors = color_list.split()
    if len(colors) < 2:
        print 'failed! please provide at least 2 colors'
        return
 
    colvec = [cmd.get_color_tuple(i) for i in colors]
    parts = len(colvec) - 1
 
    count_expr = 'index'
    expression = {'pc': 'partial_charge', 'fc': 'formal_charge',
            'count': count_expr}.get(expression, expression)
    minmax_expr = {'resi': 'resv'}.get(expression, expression)
    discrete_expr = ['index', 'resi']
 
    if cmd.count_atoms(selection) == 0:
        print 'empty selection'
        return
 
    if None in [minimum, maximum]:
        stored.e = list()
        cmd.iterate(selection, 'stored.e.append(%s)' % (minmax_expr))
        if minimum is None:
            minimum = min(stored.e)
        if maximum is None:
            maximum = max(stored.e)
    minimum, maximum = float(minimum), float(maximum)
    print ' Spectrum: range (%.5f to %.5f)' % (minimum, maximum)
 
    if maximum == minimum:
        print 'no spectrum possible, only equal %s values' % (expression)
        return
 
    if expression in discrete_expr:
        val_range = int(maximum - minimum + 1)
    else:
        val_range = maximum - minimum
        cmd.color(colors[0], selection)
 
    steps = 60 / parts
    steps_total = steps * parts
 
    val_start = minimum
    for p in range(parts):
        for i in range(steps):
            ii = float(i)/steps
            col_list = [colvec[p+1][j] * ii + colvec[p][j] * (1.0 - ii) for j in range(3)]
            col_name = '0x%02x%02x%02x' % tuple(i * 255 for i in col_list)
            val_end = val_range * (i + 1 + p * steps) / steps_total + minimum
            if expression in discrete_expr:
                cmd.color(col_name, '(%s) and %s %d-%d' % (selection, expression, val_start, val_end))
            else:
                cmd.color(col_name, '(%s) and %s > %f' % (selection, expression, val_start))
            val_start = val_end
 
cmd.extend('spectrumany', spectrumany)

