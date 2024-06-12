ccp4_contact.py:32: SyntaxWarning: invalid escape sequence '\S'
  conParser = re.compile("(\S*)\s*(\d+)([A-Z])\s*(\w+)")
ccp4_ncont.py:34: SyntaxWarning: invalid escape sequence '\s'
  conParser = re.compile("\s*/(\d+)/([a-zA-Z0-9]*)/\s*(\d+).*?/\s*([a-zA-Z0-9]*).*?:")
flatten_obj.py:280: SyntaxWarning: invalid escape sequence '\d'
  statere = re.compile("^%s_(.*)_(\d+)$" % metaprefix) # matches split object names
forster_distance_calculator.py:1: SyntaxWarning: invalid escape sequence '\,'
  '''
forster_distance_calculator.py:80: SyntaxWarning: invalid escape sequence '\P'
  print('Part of LaTeX: C:\Program Files (x86)\MiKTeX 2.9\miktex' + "\\" + "bin")
frame_slider.py:: NeedRefactor
plot_noe.py:59: SyntaxWarning: invalid escape sequence '\s'
  rgx_assi = re.compile("\s*[asignASIGN]+\s+.*")
plot_noe.py:60: SyntaxWarning: invalid escape sequence '\s'
  rgx_or = re.compile("\s*[orOR]+\s+.*")
plot_noe.py:61: SyntaxWarning: invalid escape sequence '\('
  rgx_bracket = re.compile("(\(|\))")
plot_noe.py:62: SyntaxWarning: invalid escape sequence '\w'
  rgx_tick = re.compile("(\w)\'")
plot_noe.py:71: SyntaxWarning: invalid escape sequence '\g'
  restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
plot_noe.py:76: SyntaxWarning: invalid escape sequence '\g'
  restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
plot_noe.py:81: SyntaxWarning: invalid escape sequence '\g'
  restraint_block.append(rgx_tick.sub("\g<1>^", rgx_bracket.sub(" ", restraint_line)))
propka.py:1: SyntaxWarning: invalid escape sequence '\m'
  '''
propka.py:225: SyntaxWarning: invalid escape sequence '\R'
  assert PDB not in ['NIL'], "You always have to provide PDB path. Example: PDB=.\Results_propka\4ins2011.pdb"
propka.py:520: SyntaxWarning: invalid escape sequence '\R'
  Newdir = os.getcwd() + "\Results_propka\\"
transformations.py:33: SyntaxWarning: invalid escape sequence '\*'
  """Homogeneous Transformation Matrices and Quaternions.
transformations.py:886: SyntaxWarning: invalid escape sequence '\*'
  """Return affine transform matrix to register two point sets.
transformations.py:995: SyntaxWarning: invalid escape sequence '\*'
  """Return matrix to transform given 3D point set into second point set.
wfmesh.py:45: SyntaxWarning: invalid escape sequence '\s'
  dat = re.split("\s+", line)
plugins/Caver2_1_2.py:548: SyntaxWarning: invalid escape sequence '\C'
  commandXYZ = "java %s -jar \"%s\Caver2_1.jar\" \"%s\" %f %f %f %s \"%s\" %d %d %d %s" % (JOPTS, self.pymollocation.getvalue(), input, float(self.xlocvar.get()), float(self.ylocvar.get()), float(self.zlocvar.get()), self.tunnels.getvalue(), outdir, self.varremovewater.get(), 0, self.methodvar.get(), "tun_" + self.whichModelSelect)
plugins/Caver2_1_2.py:: NeedRefactor
plugins/SuperSymPlugin.py:: NeedRefactor
plugins/annocryst.py:: NeedRefactor
plugins/apbsplugin.py:1459: SyntaxWarning: invalid escape sequence '\d'
  """
plugins/apbsplugin.py:1506: SyntaxWarning: invalid escape sequence '\d'
  unassigned = re.compile('REMARK   5 *(\d+) \w* in').findall(f.read())  # Text contains PQR output string
plugins/apbsplugin.py:: NeedRefactor
plugins/autodock_plugin.py:: NeedRefactor
plugins/bnitools.py:: NeedRefactor
plugins/castp.py:: NeedRefactor
plugins/colorama.py:: NeedRefactor
plugins/contact_map_visualizer.py:: NeedRefactor
plugins/dehydron.py:: NeedRefactor
plugins/dssp_stride.py:: NeedRefactor
plugins/emovie.py:: NeedRefactor
plugins/lisica.py:: NeedRefactor
plugins/mole.py:: NeedRefactor
plugins/msms.py:: NeedRefactor
plugins/mtsslDockGui.py:: NeedRefactor
plugins/mtsslPlotter.py:: NeedRefactor
plugins/mtsslTrilaterate.py:: NeedRefactor
plugins/mtsslWizard.py:: NeedRefactor
plugins/optimize.py:: NeedRefactor
plugins/pyanm.py:: NeedRefactor
plugins/pytms.py:: NeedRefactor
plugins/rendering_plugin.py:: NeedRefactor
plugins/resicolor_plugin.py:: NeedRefactor
plugins/show_contacts.py:: NeedRefactor