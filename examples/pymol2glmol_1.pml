reinitialize
import pymol2glmol

fetch 1TQN, async=0
preset.pretty_solv("1TQN")
select heme, organic

pymol2glmol 1TQN
import webbrowser
webbrowser.open("1TQN.html")
