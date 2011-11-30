import pymol
pymol.stored_ss = []
cmd.iterate('all', 'pymol.stored_ss.append(string.ljust(ss,1))')
print string.join(pymol.stored_ss)

