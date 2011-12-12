reinitialize
import propka
cmd.bg_color("white")
cmd.set("auto_zoom","off")

results = []
python
#resis = [["1DSB","*","30"],["1ERT","*","32"],["2TRX","*","32"],["2TRX","*","35"],["1EGO","*","11"],["1EGO","*","14"],["1MEK","*","36"],["1IUE","*","283"],
#["1PPO","*","25"],["1MEG","*","25"],["1QLP","*","232"]]
resis = [["1DSB","*","30"],["1ERT","*","32"]]

for p,c,r in resis:
	cmd.fetch(p,async="0")
	cmd.refresh()
	pkavalues = propka.propka(molecule=p,chain=c,resi=r,logtime="",makebonds="no")
	results.append(pkavalues)
	cmd.refresh()
python end

python
for p,c,r in resis:
	cmd.enable("%s"%(p))
	cmd.show_as("cartoon","%s"%(p))
	cmd.select("%s%s"%(p,r),"byres (%s and chain %s and resi %s and resn CYS)"%(p,c,r))
	cmd.show("sticks","%s%s"%(p,r))
python end
cmd.zoom("all")
