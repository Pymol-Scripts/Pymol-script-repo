# bni-tools Copyright Notice
# =====================================
#
# The bni-tools source code is copyrighted, but you can freely
# use and copy it as long as you don't change or remove any of the
# Copyright notices.  bni-tools is made available under the
# following open-source license terms:
# ----------------------------------------------------------------------
# bni-tools is Copyright (C) 2008-2016 Georg Steinkellner / Schrodinger LLC
#                         All Rights Reserved
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its built-in documentation for any
# purpose and without fee is hereby granted, provided that the above
# copyright notice appears in all copies and that both the copyright
# notice and this permission notice appear in supporting documentation,
# and that the names of the authors not be
# used in advertising or publicity pertaining to distribution of the
# software without specific, written prior permission.
# The authors DISCLAIM ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS.  IN NO EVENT SHALL THE AUTHORS BE
# LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
# for plugin description begin
'''The BNI (Beyond Normal Interaction)- Tools is a 
plug in for the PyMOL molecular visualization 
system which adds additional functionalities and 
presets to the PyMOL GUI. 
'''

from __future__ import print_function

# for plugin description end
from pymol import cmd
from pymol import util
from pymol import preset
from pymol import stored
from pymol import menu
from chempy import cpv
import os
import base64
import zlib
import sys

if sys.version_info[0] < 3:
    from Tkinter import *
    import Tkinter
    import tkSimpleDialog
    import tkFileDialog
else:
    from tkinter import *
    import tkinter as Tkinter
    from tkinter import simpledialog as tkSimpleDialog
    from tkinter import filedialog as tkFileDialog

__version__ = "0.39"
__bni_year__ = "2016"
__copyright__ = "(c) %s" % __bni_year__
__pymol_version__ = cmd.get_version()[1]
# Look for integrated BNI Tools in menu.py by checking preset settings
show_bni = 1
menu_list = [1 for i in menu.presets(cmd, "None") if "track main chain" in i]
if menu_list:
    show_bni = 0
#
# Pymol Tools v.039
#
# AUTHOR: Georg Steinkellner
#
#.
sidebar_stdwidth = int(cmd.get("internal_gui_width"))
# set keys
if __pymol_version__: #<= 1.41:
    try:
        #cmd.set_key('CTRL-X',cmd.create,(None,"sele",-1,-1,0))
        cmd.set_key('CTRL-A', cmd.select, ("sele", "visible", 1))
    except:
        pass
####################
# ICON and IMAGE GENERATION START
###################
#DATA BEGIN
bni_gif = b'''
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
MoZSzsk5vTB+GA4Ddv8CMmpssw==
'''
about_gif = b'''
eJwVlGkg04/jgHeZbXZ8FHK3SiVH30lKIUO5ksZXhNQmQsSU5IjmZq45kyMjRUZNzlzNfWtz5whR
lGREzvb9//5vnjfP6+cxMjHUPEdRARuDFv4DgdD/ie5Pwsg9RUvswY7msITUDyi9QP7zCqRRDDpb
qnJEIKHOgp0rBOlVHdCpUFYToC/WApcqquUvg0wbRcxbQNe69XX3tM1Y+GuNWNt+nF0n6CYPROmv
Oetw3LZe+mbbWfsakNsI4DwoeqdP6U4n3n2QRK6wcrnj6ESxthU433E65sW1dXe/4+rs6uZCcq2z
v3evzsJL9f6Qg5fXXQ83Kw/OLW+qh6e7gv/4KV+eObXZ8aHPHdf/nB75unjXPXh438PPz9y/xyXA
3yhwQOfJsGdAoGr49CN/X2+fWr3QMe+gYI8nwZTHrQ+ePHkcHOgVFuobEmoRMfokJDg4LDIkjHY5
diIgMjI8IsyXHhMWRQ+KodvEj3U9iI2KiYygx0fToyJjE2Jio0Pi45sD4qPiEulxMQ7JU9HxjNh4
elgiIy4hNpyRRE9MikxOSUxK8EwdTkhOs8j5lpTCiElP74xIj83ISE1PSX76LD5pJ/VZNjn/a0bm
U0ZOjjdzJiMnNysnMyWXGZA/lZmb71n4NT3/RVqmICe/IJU5yMzPzXj5MrOwMOT1fH5h8cvCgucs
VkzJ7Mvi0kfvfrxmFca9mY96t/CCzS59W5L2bq6krJz97m1RRWV2xfSb8sqiN4LkmgV2Zc10UWVl
dXlJ9TS7+lNFbX1NbVV1/Yfa+vf5Lct1nJYyzuyHpsamFg6rY+lDW0dre8uH5p3mju76ru72zraF
D/2d3R3tPf2cvv6urpnu3s7R7p6ufm4Nd6W3v6el/3sPd6Cf29s3MMwb/PhxYK91dHVwZIA38ql/
ZGlkbKhpan1sfGRkYurTxOjo5HT/1Or45KeJqfHJzxNLn2Y/Te18npmanv08+WV+aG5rdm6G/3lh
bP73ly/LQ4sb8/Mr83O789/mPi+sfFv8OvV9c+HHt68/lr8vLX77sft5ZWNp+cfyytIif3uJ//sX
f3ll7dfy2sbqb/7S+t7a+urvjbX1zd9rm9t/ttb/h/Wdva2dzT+7gu3drZ297d2/O3v/7Qr+2wOB
QIcEgOLtE/gHPs5unu7OHh6ezvdBh7ZhYBToP5DK/zzo/zsBIQSgv0hZYzazKQyxX+lGJzuvJQo4
YpptXJbfFiupfneys+xFRyLegC5r8s4AjBxY356zZIIU0xRud2eblL+6ElVhJzEzkJDbj/1O/ixr
WlHEe4Eu/16gbp4bnq7GHjYwrSweZuG3Q77oG0mdDxqoieOJ98tbv1AkZTv0nMlv8JnlXhke9Zys
8dsb/2zGqn7b5NOg76x0NC6HFDp7F4m9ppzN9eZlVE7Y8B0qXBOEHGNln0lV3o8SgZxMOzF4gjL3
zUK+862dyx12WYO9A/fh6UB/3sPnN0gHNe4Mr389umU+9nrg5X2KlseKhsednKwd94XsZx7j6O2d
jZ63YpknHRoFft8dF65pj08rpEnhGmMzLfX+eu8MHLWpa2mu9esSs7r1nP/z/fFO5TXy8asnPOOM
bvndz3YcI/7u5N0MpHQeMlGrwTglyv2iqCfq5tIWIhQ2RU+bSzlshWAvAM8PG9FmIngI9Y35a1do
BF9w7oVKUOfNiOSLYR+pKLdqBwPDkfWuA2b5eYaSlMeUmCgI8gU85Vb4rGGq7akBq2Nm/WvxtlBk
i9Oxoc6wq/LXq7oKFHA76gkHbwTntaddl17DxAo6bxtlzVY+67KB3NuKkG1jl19VvNc4MEicY65G
4HmtR0jnIyv8XfwsavX6ooZrvt9RdBuxLMy6iziCS0IctSXZCdxez1xULHnztdhsiPi6/UQDSD32
YDIi5ZUujylqKOXLR6EJL1NFDRV59SkmkuR6y74o51BRw2O99WG4HM+aL93KJcYQQyny0Kln2j6R
IYbHKlXsYwNI0ocdtP1TZw1P7KU7xsg5uU4gXAUe1yNlm5nliqkd401x2tCI3wX47vo53IVQmreh
VHqriNnJXKZ3u0pyJdxZ/tX6hJdq2siXEoWy+jXTE/7E78JnnI75LBHZf1RMbmmlpJpKuVw/jFaf
Dm8xU+OtH0dIuSIGRAan48L6TpYjZHDHXRkU1SoI+IEhbyo4MFKWo2iA0eBRT/1J3WDcbz9Adbrn
oDlr1dItaHdiJ/dfWBvyECZ0atuYlTunJtHzN+ZD4M2oJ/fhkncZq0LHBraPqEmJ2WtXpL+RdcB8
HkbcvqxmZW0Al4pH/Lv1gBt/O1HVmaEnSR0SOyIs5YL4XdR/t1K6XiePSYVym0atUJJON6ITdLY6
B/u1lITKa2Ua8Pv7tf4s3r6q5Tfdgfx9a7+TjlSGj/tV+ZucQ2aaq71HEpVZl5t7Nak0sKFSs9Wg
uc4Ox0hnTJQ+vD377yo77MAvWuaj1iPyK5LCJzIWCP5DUDJjouO432bMVBJFInEnUrqfe05XbeX0
MFkYtGGRWBs3J+cfY8hBW8TVafND1jh9KhAJlHPBbXnBb7+og9ajsb1InB+UJ4yyYzNsHofAlTRD
FfASWLEajiM0KUqEBGCfLB4Rs3TVEykAUHkFL+UDqU2YYPiJczU0eHSREdfQMXa4LESsSRShGxmG
/EvJl3Q+hNDlqyk8WRdFlzANtdebjIPGIOibxrnVayC4pticl2BG+EgnYHJnXf9JMAiuQTZmBD2G
YsFQw/Pj3lGmqCNyFweMztshI20LoLdcqLnNmqIEmwzHh77R7VGNTfD+ayWS2fVQQn6bUjCMqsnj
5h5SimTYwOLEetmpuhlNcAIbLDk7Myja4ZZoyeRgUJ8MPg58Ccd9Ags9wxuvzlPLTpbsF6IIaCZ2
89zyfxHzcrqLz2S+Ke7b03sXFMcW8oJ9Ef77Xv8UKin8gIRflPZ8Xr4/reXickeDll8IOhcfYSyk
ZvFP/WHpgozFg2qfZAzMTovPZjjv3oyRCyIeFrmJFzoznqTr8GuO8Vox8mA9CN77HSrjwQ89qJRO
j9HPl0Tt7zxJ0LtW/CvfkJehKItuin4TsCGj5C1kIXTI0pJ4+PVNZp4o7lTHw+IZ3LJxpxKJGn6z
0+Agj73wsSgk9CDH6SDYOVKXbhRmO/RCvpwJRukfYlTcyFOJJ4LFYYeRtYYGB3s1IXoG+YlkziER
ppCq7A21RJ3P4IrymnR9T+nxkap00R9ZgjB3+1gk5fGLfcw0BX1JsJ0lCSqX6wy9IhIaM3lGtYaj
2HGBmiX078UXB26CSpXHKfGL3/Ic2/pCPYWaJ0b0D8tuuLZ72O/GMmsc6/OPSkkkKdY//8vJRbyH
oPL1GD7ca/i7xCZl2CkRMR89vDmxSQUxcOvEaQmhYcBYqwRBFkO8MMvCmOjkr8ZouzahJJ/rKxWd
FUbXQVxQlFaxEm7kucnZ1uOeLdGU3YSAojEUs40tjtBDURpaRK0z9I++i2FcjNvaz5OMoo8SboTO
+KG+G0TR2xiTQQRIS6OzsAJhMrLHlYJ6wPlFr9MAjGh/uq1p97TTQyiqFbO9rp4GAtXMZBR6TqzX
nROuylDDeieklDapH9W9HMJQ/bxfSAlAnCLNwEI6D0hC3MCoV6eAc6d4B90OJfuK5agHXtRU3iAK
a3isOikSKZjczLBCjUPiEQrNgxonI09Zl6It/TKBqcwoSSY5wqpXf4hUdE/jzm/k1emQHvhshzQp
/pmtu17a2KrwG90mtI9k3QHmOQ2xC3md3r3SKDq5XYzJjaS/l0Y9AFV9Cy2LNfU7P7FCiCxcIFQ6
1kxMXRnsqEaqxXqpbwfHnbRcQycXPFoXP0QXPFj7aYu+AwTg6f3VurJB5FAzaaG0xtMuC+ALNmfz
kHWuRoczmeEaF6GOyNmWNsEzh4RHGXtyBaTQfVvcUj9XPSSJ/ExZEytkTx6QOG2zqOer93sooQVl
yoYQ8s8zbhFnv7zZKhXPJEc40CIBpY6LmetJqISzefuIBwzagvPuK7YefgPJu3zuQl74Mc6OXJPi
rvTdsYSq6il7vKX/YbpF/LkbssA0STNw1Sz76OFwlCInUp8WAtGoTxGHIqmqiDgI1CBgehzJD2Mp
rCbNvJAlgSCAowXW3lUekUB+TfxGQc9uvcBmvm0tROij+qIm5q1BofvurDafqbl5fiVaB/jGkArt
NRDJ9DLAdCVtH1sPXRTgH/DdMxLGlJskmLMGP6I0Es+53kReAoWKKwmO4WageugxYOHC+iqybzTP
vo5bPzJqzzB1ff7159ilP7HsWbvVoOw9hLCunBdgvb4h/iQMN8tG2RAB8IiNvoX2WZJIVNr12GGq
9n7ttIewp2m2aaOij/HwBhukFhtqSrq4iofLjUFjOJbzeLgOEUGuhMstIlhngUxuiDtXdYUL02L7
f6Xhvo+GZHL391RCz7Ch7my0fhnajhtCpF1ZIWFNqWGVjLBlvggNBNOhhr3mqi5zw4kEESYPLleE
OPjUhMcNl6uBXiVeySVFOpLRvukIua2AaiB8QtgXoidDB4W/jcI94IZYgxx7Gchkmt9eldfb+1eV
CP4CYVHbJ+1+OE1gi3x1aV1oB+bfBCVBHobI8OBnyNA6YcR1dswZPPxpuI0pGapOgJpy4eFsuD4Z
SiBA6RysIwF6hon4RxFIo9IdaRFEYuwbrtgWIZbAhPEI4u986epcqCNH+DE/pkITSuNf2CJidYix
SoUwHTK9cAyaRjLY4mL1CVAtjv8OHv42A0FZvJRAi2sxvnLyT5wWU+gxm36ZGV+3hX3Hp59hwi5T
45VchDaTELVjqntktBafHs7ArnBQBSHYBkI0B6qsO0MQsUD78qFKXIg5DbEurqZIhJKPAFsAyhxv
WO5omIo44ck09A6PkQ4RUgEEKaekUvRnUGlfLn1gw9Tfp5hJYK2BFBpBBZqHUNuGXXeFqj1DDVDQ
BCDFkS++x8A6slPU3aA8TdPtqFSrLgiBDyugnQcGscZI6IciHKEL8piRELQGO/UeOtSKW2bDr7xP
u/EJPpWM8pFKIZBOPAZSrnOcvxY9LSNgS0gZWGSS+uyTeKoCoIYLx8NUqJBkN+h14qlwS8IcR5hx
4tlxDkSfgWgGMMv4UC1y7CGpMC1q4gobpoWHhb3et8OOKt0WSX+H+0lDXgZBs/0hqW4onw6ZRU6m
81iIPlMojXp6jq0kYGeWmiQUxCTgn5rU4WFyNbChedwuF5spYGTVdUG0aEId/qrVhOz599DLJPQy
NUqL/PS6E2boPbT+kSba/4maKG57AxO3DXyyxfQy0ZV4dLEh5jvZZYuD3bMVWc2/O+J98DcNwz6k
BaPg5gCZg2vYYeYZEXEcthhdl35gl4+d5ub+9/puN+AklG8k0X23BC43y8X6pucPEaQwKJmQ+WhG
2V0ryyh8NxqBQivY5U+oB9fm5V+0w6yScKZhp+OOXeYx0Zsp6uvn9rFDzbdf526Q7y/yMYt87G9A
esz/ZeYf+2M2R8akLk4Omv0h3y+evbTGxz7SAHbJr7yKjf6A9GxrX3pKS+9wXx1/dH+VZAwrK/S3
eyW4bae8BdK7QZfbIxVS6IVYK+Ow8KKggPuPKnE2mUUmrYXutY8siyEUBvbaedFdGq7LF6A9yD+j
VvzSVLywtlh9NjlPmvXuPCvahFUVwKqxY72vZdVLF7PHWRzpkjOrLBe9kuaAEjWetSI7A2sUfY1S
InWIOnC+dMiudCSgdCyzdLy2dBKlnG0XrcS9JOqq3G+3fyrgzWLmmx+1b36Oq0+Pl0SMlszTS9bp
4t2Zt3DIEMTqJcppQ2exkA0a4iyg7XYI94QIPcu5vkwWvttqeO0L/DQXYkNCr9m9lYgV697BWLQe
f7AevZl6VfgMaTxJ4nCWRCXIj/0OJYpQVuIKkhIflRHsxbp/4JLEr/4TC4Q7X49+eT1uw0+XhFY0
EWI4X9J3xyVp2EIowkP4yytOwnPS5ST3N1/HqatU6Eka9mwgAP7fJCzIwpHRwunkwSQJxzrySSqw
xvfaka70UH75b2B5IQ1irmWENALQm36kHb/kSEggy4g870L7U3GPBQxykZGBNdFZNfQJ5DdOlXqW
Cdga7VPn6rv0yFa55ilLjl5jnqFVI2yBgJLB+Vq1Bfa1L0VA7JkaL3eF5EAaMs8vbsGvOlmhQkYb
qiUCqOMMIaUddogTRFEco9ZitwbSdPVBAQ1dUx6vWXdiu4I+IrJFwK6katcrlzQKROqbRSpbWeeU
uGWMMfXkmip3+7q2wJih3arnS7jRJVwNGXKMDznPRJm0wjT4kGM0GFS0vsIOLE4QnVjCpUU0HrFN
lt5++yXrfV8d22td2C2w8XsWsMTCHSCChGkgSSIOUoZEUEEwPOgF48Prng97WWBlArKk+UMQCdi7
DCg/Bg7vNRM8RFVlRYd3OUFljQIqGkEDYQjgxBncOh+MmQEJk6CwvA+EuGaV3g/j7Fbzn80msi3i
ai2ARcuZSQnNnzgvLkiMBD7MxKaSwGIEsCgeayTb/HcJtwm0y7JbxRmwqyqAVTbwQKXFzaxlMbCt
OrNtjQAR54AFmnhsGAcsQwTrPm4XEFo2SWBpdqsFCVxKAKxudJBuAJt8xCkmWBoPTaXioN5dz+q7
clS63MxwQSRIXDYsae9NQD3L9ieuHtR5px4XgIecm+mKnmx2vwG2YLeOksBHaT0pjztQem2bRJiJ
Gez5494ulV63eiGTEq9Sj4YP9kAYB0MtFu4iA8Z7zVJk+I+Z1iASOBbf0WrWBddrS+R3FBFx6zSU
LAmcOoNb6gUv7/Wvan+MyW6XJoGf43GbRPQ1fN8guuStR7Ewr+8fbcjzmVY+ByxNbm+ZBKxUoAwz
oIgEPjeDk2johiTxZBo+QmY7LNiIVvKTiOsCLtyhD5qE86ZiYWqIek6fqwqwTcVcI4FTiN2QWaQ4
gbu/D7Ca/CipDaxzOgXEIZj3ULzKAFv2414dpJ/IVVruVpkaPEyGzf7EJeEhp0hD+J9DkDysOAGi
nWMi0Sei1wAYOIDJmFbj0iHD+CHHIFwpVzhYduDRsRJrFZatClaaDbYtxSYBEBkSOIQETJkBTjlg
F0yPOJm3YNbrY9YrF8+T7uMd0AHw9/rheiOpOQOB8R2BetFBcixbWRhAAsuQID8IiCgdMIYwTJlq
3mNizpWOWu11ujv0tt74WDrTX8r4uL7X0XQF16L68fkyzoIj8uvn4BWdcXHB0gd/9L7CeJamCkpk
BrSfCMaUYmXxEwgQRD+PS9Ae7FEFvjDADQ6fYN44A27vXdlpQklvT07rHeJwrel4i92URV3jEgcE
Z08pgPY/nQHJmIGiSBMiDkKtQAdErT/lCi8V052pAk2h9dyN6zrwfNalT+Qa8PEogJCWf5mvqu06
xFZeKoXMXiQpXOL5G0JktJtBZfsul1HtUFYzEFW+ULTv1Z3ut44XDRWffXFPZS0tQcF4kNhPEZN5
MJzbBGFjIHkfWupqF+uu/vUh//qB+1Fhviz2TuWXtdxr3Hbbe2JwuXydOX7Xrz1Lgu/MWqrDTHD/
1xNWaBY7w8UIfH5hrz1vUNv9ELF13VDZNkPL5uFDG/8IiEeqhK6ymGlPxX1l4Ns8hEXU6SVefRpY
fm6Cph1IDg+uGAsaBihfx2O/7lhXoauq1r6+E5BvGfErkUlVunjAknnrmLw/RKKOwiKzTi4UCIBX
8otFWO2guirakrUw/OUev1olyzW1bhpXSH0k45azNIriQMcblyc/L09dANGw5aEXpmyzFMDV0VmN
7x1uVg/Ia/88CQx/rq86j/YjrfwJXtl6vuI3gynu6c6oq9uF+4WeHa/Oqhv3XroR+AVA1DlM+M0t
TV+7x5oTrHSpszV/4QSbUJBmMfQ6LUtQz9NqSBUkldBCv2d6rorH1j9vrF+TX20udhUI+FC1hkxl
6pfdbwsFcjskuJv5b5iMtju2McX410lciS4RqpbIgSZhfQnsrttJkq+TZE80osygaA+FW1bJJ2fm
ng6sH8nShnokw0vWL/1TImv5xxVXUs4td5v+U4264Hlws+KH2f1bm75PNv1yN7GiAEjr/wCgJh/y
'''
# DATA END
def getpic(picture, decode=0):
    if decode == 1:
        pic = zlib.decompress(base64.decodestring(picture))
    else:
        pic = base64.encodestring(
            zlib.decompress(base64.decodestring(picture)))
    return pic
bni_image = getpic(bni_gif)
try:
    bni_image = Tkinter.PhotoImage(data=bni_image)
except:
    bni_image = ""
def do_nothing():
    pass
####################
# ICON and IMAGE GENERATION END
###################
###################
# ABOUT SECTION BEGIN
###################
######################################################
bni_info = '''The BNI (Beyond Normal Interaction)- Tools is a 
plug in for the PyMOL molecular visualization 
system which adds additional functionalities and 
presets to the PyMOL GUI. 
Version: %s (%s)
Please cite BNI-Tools if you find the following 
additional features useful. 
   # BNI tools dropdown menu 
     [+]Load Files->
      [-]modeller files
        * load multiple modeller files and sort by 
        * objective function energy 
        * and add energy to object title
      [-]autodock files
        * load autodock .dlg or .dlg.pdb file into
        * different states with energy in object 
        * title
      [-]amber minimization
        * load amber minimization file with energy in
        * title if the .info file is present in same
        * directory and has the same name as the pdb file.
      [-]delphi phi,dx map
        * load delphi map and corresponding pdb file
        * simultaneously and show surface colored by 
        * PHI or DX map. (show the surface to see the 
        *  effect)
      [-]casox map
        * load casox map (cavity calculation) and 
        * show ligsite accessibility 
        * value maps. (7 closed cavity to 1 open)
      [-]multiple files into states
        * load multiple pdb files (e.g. MD simulation 
        * snapshots) into one state. (object is named
        * by first object loaded)
     [+]Fetch->
      [-]Density View (EDS)
        * load density and pdb file from EDS
        * (Electron Density Server)
        * if available, and show density with density
        * wizard
      [-]RCSB Biol. Assembly
        * load biological assembly from RCSB
        * protein database
      [-]2FoFc map(s)
        * load (multiple) 2FoFc maps
        * from EDS density server
        * if available
      [-]FoFc maps(s)
        * load (multiple) FoFc maps
        * from EDS density server
        * if available
     [+]Edit->
      [-]HIS --> HID,HIE,HIP
        * change histidine residues to HID,HIE,HIP
        * depending on hydrogens on histidine
      [-]HID,HIE,HIP --> HIS
        * change altered histidine residues 
        * back to HIS
      [-]Poly-Alanine Chain
        * create a poly alanine chain (GLY and ALA)
        * for molecular replacement
      [-]MSE --> MET
        * change selenomethionine to methionine
      [-]del alternates
        * delete alternates in selection
      [-]Unbond- >
        * unbond atoms in selection
     [+]Images->
       * create ray traced images depending on
       * x size and resolution (dpi)
     [+]Create->
       * create compiled graphics objects (CGO)
       * these objects can be altered in color or
       * transparency, and they can be dragged
       * and rotated in space by the
       * "action->drag" command and using
       * "shift" and mouse buttons.
      [-]Plane
        * create a plane (with certain cushion)
        * using a selection of three atoms
      [-]Box
        * create a box around selection
        * the whole box can be altered as group
        * or by side planes separately
      [-]Triangle
        * create a triangle using a selection of three
        * atoms
     [+]pseudo center atom
         * create pseudo atom at the center of atoms
         * in selection
         * move atoms with editing->"shift"
         * and middle mouse button
   # BNI tools integrated in PyMOL sidebar
     [+]Action on "all"
      [-]delete enabled
        * delete all enabled objects or 
        * selections
      [-]invert enabled/disabled
        * disable currently enabled objects
      [-]combine selections
        * combine all enabled or disabled 
        * selections to selection (sele)
     [+]action
      [-]sequence
        * show sequence in different formats
        * of selection
        * copy and paste to text file for later use
        [.]fasta
          * show sequence of selection
          * in fasta format
        [.]pir
          * show sequence of selection
          * in pir format
        [.]modeller
          * show sequence of selection
          * in modeller pir format
        [.]list
          * create residue or atom lists
          * of selection
     [+]preset
      [-]track main chain
        * create a new object which tracks the
        * main chain atoms and shows main chain
        * and side chain polar contacts
      [-]symmetry surface
        * create a symmetry view of the selected
        * atoms showing the contact surface as well
        * as a selection entry of the atoms in 
        * contact with symmetry mates
        * (only includes atoms of the initial selection).
      [-]hydrophobic residues
        * show hydrophobic residues
        * depending on hydrophobic residue scales
        * by KandD (Kyte & Doolittle 
        *           J Mol Biol 157:105, 1982)
        *    Rose  (Rose et al. 
        *           Science, 229, 834-838,1985)
        *    GES   (Engelman Engelman et al. 
        *           Annu Rev Biophys Biophys Chem,
        *           15, 321-353(1986)
        * (no window selection; just raw categories
        * are colored by: blue-hydrophile
        *                green-neutral
        *                red-hydrophobe )
      [-]surface inspection
        * create selections for surface
        * inspections
   # BNI tools additional settings in PyMOL sidebar
      [+]color
       [-]by ss
         * color helix,sheet and loop separately
         [.]Helix
         [.]Sheet
         [.]Loop
       [-]surface
         * color surface separately from atoms
         [.]by atom
           * set surface color to standard
         [.]by map
           * if a map or ramp is loaded
           * color surface by ramp/map
         [.]my color
           * color surface by own defined color
       [-]mesh
         * color mesh separately from atoms
         [.]by atom
           * set mesh color to standard
         [.]by map
           * if a map or ramp is loaded
           * color mesh by ramp/map
         [.]my color
           * color mesh by own defined color
       [-]label
         * color labels separately
         [.]by atom
           * set label color by atom
         [.]my color
           * color labels by own defined color
       [-]stick
         * color sticks separately
         [.]standard
           * set stick color to standard
         [.]my color
           * color sticks by own defined color
       [-]my colors
         * use/append own defined colors
         * own colors can be defined by
         * Setting->Colors..->New
         * to keep color settings for
         * other pymol sessions
         * you have to set the colors
         * in .pymolrc or similar
         * pymol setting file
         * like
         * set_color mycolor,[ 1.00, 1.00, 1.00]
      [+]show/hide
       [-]surface flag
         * set surface flag of atoms to show hide
         * or ignore fro surface calculation
       [-]transparency
         * set different transparency types
         * on selections or atoms
''' % (__version__, __bni_year__)
class AboutGui:
    def __init__(self, parent):
        self.buttons = ['Ok', 'Cite', 'Thanks', 'Info', 'Exit']
        self.lastinterior = None
        self.txtcolors = (("Green", ("#")),
                          ("DarkBlue", ("[+]",
                                        )),
                          ("Blue", ("[-]", "[.]"
                                    )),
                          ("Purple", ("*",
                                      )),
                          ("Red", ("ACKNOWLEDGEMENTS",
                                   )),
                          )
        about_image = getpic(about_gif)
        try:
            self.about_image = Tkinter.PhotoImage(data=about_image)
        except:
            self.about_image = ""
        self.dialog = Pmw.Dialog(parent,
                                 title='BNI-Tools (c) SCHRODINGER',
                                 buttons=self.buttons,
                                 defaultbutton=self.buttons[0],
                                 command=self.execute,
                                 )
        self.defaultwindow()
    def write_info(self, bni_info):
        listone = bni_info.split("\n")
        if self.txt:
            for line in listone:
                style = self.check_syntax(str(line), use=1, case=1)
                self.txt.insert('end', str(line) + "\n", style)
    def get_info(self, bni_info, joinit=0):
        a = []
        s = ""
        if joinit:
            a.extend(b.strip() for b in bni_info.split("\n"))
            if joinit == 1:
                sign = " "
            else:
                sign = "\n"
            s = sign.join(a)
        if joinit == 0:
            return "%s" % bni_info
        else:
            return "%s" % s
    def check_syntax(self, line, use=1, case=1):
        '''Check syntax highlighting.
        line: str:line to check for syntax
        use:0: just return empty string.
        case:1: case sensitive
        '''
        syntax = [""]
        if use == 0:
            return tuple(syntax)
        if case == 1:
            low = str.lower 
        else:
            low = lambda x: x
        for color, tags in self.txtcolors:
            for tag in tags:
                if low(tag) in low(line):
                    syntax.append(color)
        return tuple(syntax)
    def tag_list(self, txt_parent):
        '''Tag list for tags to color log output.
        '''
        for color, filter in self.txtcolors:
            txt_parent.tag_configure(color, foreground=color)
    def destroylast(self):
        if self.lastinterior is not None:
            self.lastinterior.destroy()
    def defaultwindow(self):
        self.destroylast()
        self.w = Tkinter.Label(self.dialog.interior(),
                               text='Beyond Normal Interaction - Tools',
                               background='white',
                               foreground='black',
                               image=self.about_image,
                               pady=20)
        self.w.pack(expand=1, fill='both', padx=4, pady=4)
        self.lastinterior = self.w
    def message(self, parent, text="", hull_width=400, hull_height=90):
        fixedFont = Pmw.logicalfont('Fixed')
        self.txt = Pmw.ScrolledText(parent,
                                    labelpos=None,
                                    columnheader=0,
                                    rowheader=0,
                                    rowcolumnheader=0,
                                    usehullsize=1,
                                    hull_width=hull_width,
                                    hull_height=hull_height,
                                    text_wrap='word',
                                    text_font=fixedFont,
                                    text_padx=4,
                                    text_pady=4,
                                    )
        self.txt.pack(padx=0, pady=0, fill='both', expand=1)
        self.tag_list(self.txt)
        self.txt.appendtext(text)
    def showcite(self, text="", label_text="", hull_width=400, hull_height=300, vertflex="expand", horizflex="expand"):
        self.destroylast()
        self.sf = Pmw.ScrolledFrame(self.dialog.interior(),
                                    labelpos='n',
                                    label_text=label_text,
                                    usehullsize=1,
                                    hull_width=hull_width,
                                    hull_height=hull_height,
                                    vertflex=vertflex,
                                    horizflex=horizflex,
                                    )
        pixelperchar = 12
        hull_width = int(len(text) / 3) * pixelperchar
        hull_height = int(3 * pixelperchar)
        self.message(
            self.sf.interior(), text=text, hull_width=hull_width, hull_height=hull_height)
        self.sf.pack(padx=5, pady=5, fill='both', expand=1)
        self.lastinterior = self.sf
    def execute(self, action):
        if action not in self.buttons:
            self.dialog.destroy()
        elif action == self.buttons[0]:
            self.defaultwindow()
        elif action == self.buttons[-1]:
            self.dialog.destroy()
        elif action == self.buttons[1]:
            citetext = "Steinkellner, G.(%s). BNI-Tools (Version %s)[Computer software]. Schrodinger, LLC." % (
                __bni_year__, __version__)
            self.showcite(label_text="How to cite BNI - Tools...",
                          hull_width=400,
                          hull_height=120,
                          horizflex="elastic",
                          text=citetext)
            print(citetext)
        elif action == self.buttons[2]:
            namelist = ["Schrodinger LLC",
                        "J. Vertrees",
                        "T. Holder",
                        "K. Gruber",
                        "M. Uhl",
                        "G. Oberdorfer",
                        "A. Lyskowski",
                        "C. C. Gruber",
                        "the strubi group"]
            text = "Thanks to ...\n%s" % "\n".join(namelist)
            self.showcite(label_text="Thanks for contributing, support and testing...",
                          hull_width=400,
                          hull_height=120,
                          horizflex="elastic",
                          text=text)
            #print(text)
        elif action == self.buttons[3]:
            self.showcite(label_text="Info about BNI-Tools",
                          hull_width=570,
                          hull_height=120,
                          horizflex="elastic",
                          text='')
            self.write_info(bni_info)
        #pass
###################
# ABOUT SECTION END
###################
def bni_toggleHelp(a):
    print('Toggle value:', a.bni_toggleHelp.get())
####################
# system variable settings
####################
def get_env(variable="BNI_SHOW"):
    '''Get the environment variable if set.
       variable: variable name
    '''
    try:
        value = os.environ.get(variable)
    except:
        value=""
    return value
show_bni_env = get_env("BNI_SHOW")
if show_bni_env and show_bni == 0:
    show_bni = show_bni_env
def __init__(self):
    self.__class__._gui = self
    seletext = "act on selection (sele)"
    self.menuBar.addcascademenu('Plugin', 'bni-tools', 'BNI Tools',
                                label='BNI Tools', tearoff=True)
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    self.menuBar.addmenuitem('bni-tools',
                             'command',
                             'BNIInfo',
                             'BNI info',
                             background="lightblue",
                             activebackground="lightblue",
                             foreground="darkblue",
                             activeforeground="darkblue",
                             label='BNI Tools v. %s %s' % (
                                 __version__, __copyright__),
                             image=bni_image,
                             command=do_nothing)
    self.menuBar.addmenuitem('bni-tools',
                             'command',
                             'BNIInfo2',
                             'BNI info2',
                             background="lightblue",
                             activebackground="lightblue",
                             foreground="blue",
                             activeforeground="blue",
                             label='v. %s %s' % (__version__, __copyright__),
                             command=do_nothing)
    self.bni_toggleHelp = Tkinter.IntVar(self.root)
    # Initialise the checkbutton to 0:
    self.bni_toggleHelp.set(0)
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    self.menuBar.addmenuitem('bni-tools',
                             'command',
                             'About',
                             'AboutHelpCite',
                             background="lightblue",
                             activebackground="blue",
                             foreground="darkblue",
                             activeforeground="white",
                             label='About / How to cite',
                             #image=bni_image,
                             command=lambda s=self: AboutGui(s.root))
    # Pseudo menue for clicking if clicking stucks
    # Symmetry BEGIN
    #import info
    if show_bni:
        self.menuBar.addmenuitem('bni-tools', 'separator', '')
        self.menuBar.addcascademenu('bni-tools', 'symmetry', statusHelp='Symmetry settings',
                                    traverseSpec=None,
                                    label='Symmetry ->'
                                    )
        self.menuBar.addmenuitem('symmetry', 'command', background="lightblue", activebackground="lightblue",
                                 foreground="blue", activeforeground="blue", label=seletext, command=do_nothing)
        self.menuBar.addmenuitem('symmetry', 'separator', '')
        self.menuBar.addcascademenu('symmetry', 'symview', 'Create Symmetry',
                                    label='Create Symmetry View ->'
                                    )
        self.distances = ["set", 3.0, 3.5, 4.0, 4.5, 5.0, 10.0]
        for item in self.distances:
            self.menuBar.addmenuitem('symview', 'command',
                                     'Distance %s' % item,
                                     label='%s' % item,
                                     command=lambda a=item, b=self: fast_sym(cutoff=a, app=b, selection="sele"))
        self.menuBar.addmenuitem('symview', 'separator', '')
        self.menuBar.addcascademenu(
            'symview', 'colorsym', 'Color surface', label="Surface Patch Color ->")
        self.menuBar.addmenuitem('colorsym', 'command',
                                 'Set Surface to Symetry representation',
                                 label='by symetry mates',
                                 command=lambda a="symmetry": set_surface(mode=a))
        self.menuBar.addmenuitem('colorsym', 'command',
                                 'Set Surface to standard representation',
                                 label='standard coloring',
                                 command=lambda a="standard": set_surface(mode=a))
    # Symmetry END
    # Settings BEGIN
    # PyMOL integration into Setting menu
    # Side Bar settings
    self.menuBar.addcascademenu('Setting', 'SideBar',
                                'SideBar settings',
                                label='Side Bar ->', tearoff=True)
    self.menuBar.addcascademenu('SideBar', 'width',
                                'SideBar settings ->',
                                label='width', tearoff=True)
    self.sidebar_lenght = ["set", 100, 200, 250, 350, 450]
    for item in self.sidebar_lenght:
        self.menuBar.addmenuitem('width', 'command',
                                 'Width %s' % item,
                                 label='%s' % item,
                                 command=lambda a=item, b=self: runsidebar(set="width", value=a, app=b))
    self.menuBar.addmenuitem('SideBar', 'command',
                             'fit sidebar to names',
                             label='fit to name',
                             command=lambda a="fit1": sidebar(set="width", value=a))
    self.menuBar.addmenuitem('SideBar', 'command',
                             'fit sidebar to titles',
                             label='fit to title',
                             command=lambda a="fit2": sidebar(set="width", value=a))
    self.menuBar.addmenuitem('SideBar', 'command',
                             'sidebar off',
                             label='off',
                             command=lambda a="fit2": sidebar(set="off", value=a))
    self.menuBar.addmenuitem('SideBar', 'command',
                             'sidebar on',
                             label='on',
                             command=lambda a="fit2": sidebar(set="on", value=a))
    self.menuBar.addcascademenu('SideBar', 'overlay',
                                'Overlay settings ->',
                                label='overlay')
    self.menuBar.addmenuitem('overlay', 'command',
                             'set  sidebar overlay on',
                             label='on',
                             command=lambda a="": sidebar(set="overlayon", value=a))
    self.menuBar.addmenuitem('overlay', 'command',
                             'set sidebar overlay off',
                             label='off',
                             command=lambda a="": sidebar(set="overlayoff", value=a))
    # Settings END
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    # Presets BEGIN
    if show_bni:
        self.menuBar.addcascademenu('bni-tools', 'presets', 'Preset Views',
                                    label='Preset Views ->', tearoff=True
                                    )
        self.menuBar.addmenuitem('presets', 'command', background="lightblue", activebackground="lightblue",
                                 foreground="blue", activeforeground="blue", label=seletext, command=do_nothing)
        self.menuBar.addmenuitem('presets', 'separator', '')
        self.menuBar.addmenuitem('presets', 'command',
                                 'Track Mainchain',
                                 label='Track MainChain',
                                 command=lambda a="track_mainchain": new_views(preset=a))
        self.menuBar.addcascademenu('presets', 'HPvalue',
                                    'Hydrophobic Values',
                                    label='Hydrophobic Values ->'
                                    )
        #
        self.HPvalues = ["KandD", "Rose", "GES"]
        for item in self.HPvalues:
            self.menuBar.addmenuitem('HPvalue', 'command',
                                     'HPvalues %s' % item,
                                     label='%s' % item,
                                     command=lambda a=item: select_hp(hp_value=a))
        #
        self.menuBar.addmenuitem('presets', 'command',
                                 'Surface Atoms',
                                 label='Atoms for Surface Inspection',
                                 command=lambda a="surf_rel_res": new_views(preset=a))
        # Presets END
        self.menuBar.addmenuitem('bni-tools', 'separator', '')
    # Load Begin
    self.menuBar.addcascademenu('bni-tools', 'Load Files', 'Load Files',
                                label='Load Files ->',
                                tearoff=True)
    if show_bni:
        self.menuBar.addmenuitem('Load Files', 'command',
                                 'Multiple Files',
                                 label='Multiple Files',
                                 command=lambda b=self: multi_pdb(b).fetchFileName())
    self.menuBar.addmenuitem('Load Files', 'command',
                             'Modeller Files',
                             label='Modeller Files',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="modeller"))
    self.menuBar.addmenuitem('Load Files', 'command',
                             'AutoDock Files',
                             label='AutoDock Files',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="autodock"))
    self.menuBar.addmenuitem('Load Files', 'command',
                             'Amber mimimization',
                             label='Amber minimization',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="ambermin"))
    self.menuBar.addmenuitem('Load Files', 'command',
                             'Delphi Map Files',
                             label='Delphi PHI,DX Map',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="delphi"))
    self.menuBar.addmenuitem('Load Files', 'command',
                             'Casox Map',
                             label='Casox Map',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="casox"))
    self.menuBar.addmenuitem('Load Files', 'command',
                             'One Object',
                             label='Multiple Files into states',
                             command=lambda b=self: multi_pdb(b).fetchFileName(ext="oneobject"))
    self.menuBar.addcascademenu('bni-tools',
                                'FetchPdb',
                                label='Fetch ->',
                                tearoff=True)
    fetch_types = ["pdb", "pdb 2fofc fofc", "pdb1", "2fofc", "fofc"]
    fetch_types_names = [
        "PDB Codes(s)", "Density View (EDS)", "RCSB Biol. Assembly", "2FoFc map(s)", "FoFc map(s)"]
    for typ, typn in zip(fetch_types, fetch_types_names):
        self.menuBar.addmenuitem('FetchPdb',
                                 'command',
                                 'Fetch %s' % (typn.strip()),
                                 label='%s' % (typn.strip()),
                                 command=lambda b=self, a=typ: bni_fetch(parent=b, typ=a))
    #Load End
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    # Edit BEGIN
    self.menuBar.addcascademenu('bni-tools',
                                'EditAtom',
                                'Edit Atoms',
                                label='Edit ->',
                                tearoff=True)
    self.menuBar.addmenuitem('EditAtom', 'command', background="lightblue", activebackground="lightblue",
                             foreground="blue", activeforeground="blue", label=seletext, command=do_nothing)
    self.menuBar.addmenuitem('EditAtom', 'separator', '')
    self.menuBar.addmenuitem('EditAtom',
                             'command',
                             'His2',
                             label="HIS --> HID,HIE,HIP",
                             command=lambda b="his2hid": his_handle(b))
    self.menuBar.addmenuitem('EditAtom',
                             'command',
                             'Hid2',
                             label="HID,HIE,HIP --> HIS",
                             command=lambda b="hid2his": his_handle(b))
    self.menuBar.addmenuitem('EditAtom',
                             'command',
                             'polyalanine',
                             label="Poly-Alanine Chain",
                             command=lambda: polyala("sele", name="polyala"))
    self.menuBar.addmenuitem('EditAtom',
                             'command',
                             'MSE->MET',
                             label="MSE-->MET",
                             command=lambda: altermse("sele"))
    self.menuBar.addmenuitem('EditAtom',
                             'command',
                             'delalt',
                             label="del alternates",
                             command=lambda: delalternate("sele"))
    self.menuBar.addcascademenu('EditAtom',
                                'UnbondAtoms',
                                'Unbond Atoms',
                                label='Unbond ->',
                                tearoff=True)
    self.menuBar.addmenuitem('UnbondAtoms',
                             'command',
                             'unbondatom',
                             label="atoms",
                             command=lambda: unbondatom(selection1="sele",
                                                        selection2="all",
                                                        selection="sele"))
    self.menuBar.addmenuitem('UnbondAtoms',
                             'command',
                             'unbondallmetals',
                             label="metals",
                             command=lambda: unbondatom(selection1="metals",
                                                        selection2="all",
                                                        selection="sele"))
    self.menuBar.addmenuitem('UnbondAtoms',
                             'command',
                             'unbondallhetatoms',
                             label="het atoms",
                             command=lambda: unbondatom(selection1="hetatm and not r. MES and not solvent",
                                                        selection2="!hetatm",
                                                        selection="sele"))
    # Edit END
    # Create BEGIN
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    self.menuBar.addcascademenu('bni-tools',
                                'Image',
                                'Create Image',
                                label='Images ->',
                                tearoff=True)
    #Ray BEGIN
    self.ray_settings = [
        "set", "standard", "journal single", "journal double", "poster", "hires"]
    ray_settings_menu = ["set (x mm/cm/inch y dpi)", "standard (10 cm 200 dpi) ", "single column (86 mm 400 dpi)",
                         "double column (178 mm 400 dpi)", "poster (30 cm 300 dpi)", "hires (12 cm 1200 dpi)"]
    for item, item_names in zip(self.ray_settings, ray_settings_menu):
        self.menuBar.addmenuitem('Image', 'command',
                                 'Ray %s' % item,
                                 label='%s' % item_names,
                                 command=lambda a=item, b=self: _set_bni_ray(setting=a, app=b))#fast_sym(cutoff=item))
    #self.menuBar.addmenuitem('Ray',
    #                         'command',
    #                         'His2',
    #                         label="HIS --> HID,HIE,HIP",
    #                         command= lambda b="his2hid": his_handle(b))
    #Ray End
    #Plane begin
    self.menuBar.addcascademenu('bni-tools',
                                'Create',
                                'Create',
                                label='Create ->',
                                tearoff=True)
    self.menuBar.addmenuitem('Create', 'command', background="lightblue", activebackground="lightblue",
                             foreground="blue", activeforeground="blue", label=seletext, command=do_nothing)
    self.menuBar.addmenuitem('Create', 'separator', '')
    self.menuBar.addcascademenu(
        'Create', 'GetPlane', 'Get Plane', label='Plane ->', tearoff=True)
    triangle_txt = "3 atoms in (sele)"
    self.menuBar.addmenuitem('GetPlane', 'command', background="lightblue", activebackground="lightblue",
                             foreground="blue", activeforeground="blue", label=triangle_txt, command=do_nothing)
    self.menuBar.addmenuitem('GetPlane', 'separator', '')
    self.plane_settings = ["set", "small", "medium", "large"]
    self.plane_settings_names = [
        "set plane size", "small (size 5.0 %s)" % "Ang", "medium (size 10.0 %s)" % "Ang", "large (size 20.0 %s)" % "Ang"]
    for item, itemnames in zip(self.plane_settings, self.plane_settings_names):
        self.menuBar.addmenuitem('GetPlane', 'command',
                                 'Plane %s' % item,
                                 label='%s' % itemnames,
                                 command=lambda a=item, b=self: get_plane_handler(action=a, app=b))
    #Plane end
    #Box begin
    self.menuBar.addcascademenu('Create', 'GetBox', 'Get Box',
                                label='Box ->', tearoff=True
                                )
    self.box_settings = ["set", "small", "medium", "large"]
    self.box_settings_names = [
        "set box cushion", "small (cushion 1.0 %s)" % "Ang", "medium (cushion 2.0 %s)" % "Ang", "large (cushion 3.0 %s)" % "Ang"]
    for item, itemname in zip(self.box_settings, self.box_settings_names):
        self.menuBar.addmenuitem('GetBox', 'command',
                                 'Box %s' % item,
                                 label='%s' % itemname,
                                 command=lambda a=item, b=self: get_box_handler(action=a, app=b))
    #Box End
    #Triangle BEGIN
    self.menuBar.addcascademenu('Create', 'Triangle',
                                'Triangle',
                                label='Triangle ->',
                                tearoff=True)
    triangle_txt = "3 atoms in (sele)"
    self.menuBar.addmenuitem('Triangle', 'command', background="lightblue", activebackground="lightblue",
                             foreground="blue", activeforeground="blue", label=triangle_txt, command=do_nothing)
    self.menuBar.addmenuitem('Triangle', 'separator', '')
    self.menuBar.addmenuitem('Triangle', 'command',
                             'create',
                             label='create ..',
                             command=lambda a="sele": get_triangle(selection=a, invert=0))
    #Triangle END
    if show_bni:
        #Set Surface BEGIN
        self.menuBar.addcascademenu('Create', 'SetSurface', 'Set surface',
                                    label='Surface ->', tearoff=True
                                    )
        self.surf_settings = [
            "Methionine", "Include", "Exclude", "Hide", "Show"]
        for item in self.surf_settings:
            self.menuBar.addmenuitem('SetSurface', 'command',
                                     'Surface %s' % item,
                                     label='%s' % item,
                                     command=lambda a=item, b=self: get_surface_handler(action=a, app=b))
        #Set Surface End
        # Create END
    # Show selection  on BNI Tools or in PyMOL directly
    if show_bni:
        selemenu = "bni-tools"
        self.menuBar.addmenuitem(selemenu, 'separator', '')
        # Delete BEGIN
        self.menuBar.addcascademenu(selemenu, 'delete',
                                    'Delete all',
                                    label='Delete enabled',
                                    tearoff=True)
        self.menuBar.addmenuitem('delete', 'command',
                                 'enabled',
                                 label='all',
                                 command=lambda: del_enabled('all'))
        self.menuBar.addmenuitem('delete', 'command',
                                 'objects',
                                 label='objects',
                                 command=lambda: del_enabled('obj'))
        self.menuBar.addmenuitem('delete', 'command',
                                 'selections',
                                 label='selections',
                                 command=lambda: del_enabled('selections'))
        self.menuBar.addcascademenu(selemenu, 'invert',
                                    'invert all',
                                    label='Invert enabled/disabled',
                                    tearoff=True)
        self.menuBar.addmenuitem('invert', 'command',
                                 'all',
                                 label='all',
                                 command=lambda: inv_enabled('all'))
        self.menuBar.addmenuitem('invert', 'command',
                                 'objects',
                                 label='objects',
                                 command=lambda: inv_enabled('obj'))
        self.menuBar.addmenuitem('invert', 'command',
                                 'selections',
                                 label='selection',
                                 command=lambda: inv_enabled('selections'))
        #combine  BEGIN
        self.menuBar.addcascademenu(selemenu, 'Combine', 'Combine selections',
                                    label='Combine sele..', tearoff=True
                                    )
        self.combine_sele = ["all", "enabled", "disabled"]
        for item in self.combine_sele:
            self.menuBar.addmenuitem('Combine', 'command',
                                     'combine %s' % item,
                                     label='%s' % item,
                                     command=lambda a=item: combine_sele(check=a))
        #combine END
        # Delete END
    self.menuBar.addmenuitem('bni-tools', 'separator', '')
    # Create Pseudo Atom for centers
    self.menuBar.addmenuitem('bni-tools', 'command',
                             'Create PseudoCenter Atom',
                             label='PseudoCenter Atom',
                             command=lambda: get_center_atom(selection='(sele)', enabled="1", onlysele="1"))
    # Get sequence
    if show_bni:
        self.menuBar.addmenuitem('bni-tools', 'separator', '')
        self.menuBar.addcascademenu('bni-tools', 'Sequences',
                                    'get sequence',
                                    label='Sequences',
                                    tearoff=True)
        self.menuBar.addmenuitem('Sequences', 'command', background="lightblue", activebackground="lightblue",
                                 foreground="blue", activeforeground="blue", label=seletext, command=do_nothing)
        self.menuBar.addmenuitem('Sequences', 'separator', '')
        self.menuBar.addmenuitem('Sequences', 'command', 'Fasta', label='fasta', command=lambda: get_sequence(
            selection='(sele)', format="fasta"))
        self.menuBar.addmenuitem('Sequences', 'command',
                                 'pir',
                                 label='pir',
                                 command=lambda: get_sequence(selection='(sele)', format="pir"))
        self.menuBar.addmenuitem('Sequences', 'command',
                                 'modeller',
                                 label='modeller',
                                 command=lambda: get_sequence(selection='(sele)', format="modeller"))
        self.menuBar.addcascademenu('Sequences', 'Normal',
                                    'get normal sequence',
                                    label='List ->',
                                    tearoff=True)
        self.menuBar.addmenuitem('Normal', 'command',
                                 'Residue',
                                 label='residues',
                                 command=lambda: get_sequence(selection='(sele)', format="normal"))
        self.menuBar.addcascademenu(
            'Normal', 'bfac', 'bfac per res', label='bfac per res ->')
        calc_fac_name = ['all', 'mainchain', 'sidechain']
        calc_fac = ['all', 'main', 'side']
        for item, itemname in zip(calc_fac, calc_fac_name):
            self.menuBar.addmenuitem('bfac', 'command',
                                     'bfac %s atoms' % itemname,
                                     label='%s atoms' % itemname,
                                     command=lambda a=item: get_sequence(selection='(sele)', format="normalb_%s" % a))
        self.menuBar.addcascademenu(
            'Normal', 'occ', 'occ per res', label='occ per res ->')
        for item, itemname in zip(calc_fac, calc_fac_name):
            self.menuBar.addmenuitem('occ', 'command',
                                     'occ %s atoms' % itemname,
                                     label='%s atoms' % itemname,
                                     command=lambda a=item: get_sequence(selection='(sele)', format="normalq_%s" % a))
        self.menuBar.addmenuitem('Normal', 'command',
                                 'Residue',
                                 label='all atoms',
                                 command=lambda: get_sequence(selection='(sele)', selecting='atom', format="normal"))
        # Print PDB line
        self.menuBar.addmenuitem("Sequences", 'command',
                                 'Get PDB lines',
                                 label='PDB Format',
                                 command=lambda: get_pdb_lines(selection='(sele)', state="0", enabled="1"))
    cmd.extend('bni_fast_sym', fast_sym)
    cmd.extend('bni_del_enabled', del_enabled)
    cmd.extend('bni_del_sym', del_sym)
    cmd.extend('bni_new_views', new_views)
    cmd.extend('bni_new_name', rename_object)
    cmd.extend('bni_get_pdb_lines', get_pdb_lines)
    cmd.extend('bni_mse', mse)
    cmd.extend('bni_get_box', get_box)
    cmd.extend('bni_get_axes', get_axes)
    cmd.extend('bni_get_point', get_point)
    cmd.extend('bni_get_atom', get_atom)
    cmd.extend('bni_get_center', getCenter)
    cmd.extend('bni_get_center_atom', get_center_atom)
    cmd.extend('bni_check_names', check_names)
    cmd.extend('bni_sidebar', sidebar)
    cmd.extend("bni_get_plane", get_plane)
    cmd.extend("bni_get_polyala", polyala)
    cmd.extend("bni_his_handle", his_handle)
    cmd.extend("bni_get_sequence", get_sequence)
    cmd.extend("bni_ray", bni_ray)
    cmd.extend("bni_combine_sele", combine_sele)
    cmd.extend("bni_unbondatom", unbondatom)
def bni_fetch(parent, typ="pdb"):
    '''Fetch PDB codes or maps.
    USAGE: bni_fetch(window,typ='pdb')
          window: parent
             typ: 'pdb'
                  'pdb 2fofc fofc'
                  'pdb1'
                  '2fofc'
                  'fofc'
    '''
    print("# -------------------#")
    print("# BNI-Tools v. %s" % __version__)
    print("# -------------------#")
    PDBcode = tkSimpleDialog.askstring('BNI Fetch',
                                       'Please enter 4-digit\npdb codes (space separated):',
                                       parent=parent.root)
    if PDBcode:
        PDBcode = PDBcode.lower()
        PDBcode = PDBcode.strip().split()
        failed = 0
        typ = typ.lower()
        if typ == "pdb 2fofc fofc":
            if len(PDBcode) > 1:
                code = PDBcode[0]
                print(
                    "# Density Wizard only for one map at once. Multiple PDB codes are given.")
                print("# First PDB code (%s) will be used." % code)
            else:
                code = PDBcode[0]
            typ = typ.strip().split()
            names_uplow = cmd.get_names()
            names = [nam.lower() for nam in names_uplow]
            if code in names:
                code_index = names.index(code)
                code = names_uplow[code_index]
                typ = ["2fofc", "fofc"]
            for source in typ:
                print("# Fetch %s of %s" % (source, code))
                cmd.fetch(code, type=source)
            if __pymol_version__ < 1.4:
                print("# Pymol < 1.4 detected.(v. %s)" % __pymol_version__)
                print("# Maps are mixed up. FoFc is 2FoFc and FoFc is 2FoFc.")
                map1 = ["w1_", 0, "%s_%s" %
                        (code, "fofc"), 1.0, "slate"]
                map2 = ["w2_", 1, "%s_%s" %
                        (code, "2fofc"), -3.0, "tv_red"]
                map3 = ["w3_", 2, "%s_%s" %
                        (code, "2fofc"), 3.0, "tv_green"]
            else:
                map1 = ["w1_", 0, "%s_%s" %
                        (code, "2fofc"), 1.0, "slate"]
                map2 = ["w2_", 1, "%s_%s" %
                        (code, "fofc"), -3.0, "tv_red"]
                map3 = ["w3_", 2, "%s_%s" %
                        (code, "fofc"), 3.0, "tv_green"]
            names = [nam.lower() for nam in cmd.get_names()]
            if map1[2].lower() not in names or map2[2].lower() not in names:
                print(map1[2], map2[2], names)
                print(
                    "# Some maps could not be loaded. Density Wizard closed.")
                return
            try:
                cmd.wizard("density")
                cmd.refresh_wizard()
                cmd.delete(map1[0] + map1[2])
                cmd.delete(map2[0] + map2[2])
                cmd.delete(map3[0] + map3[2])
                cmd.get_wizard().set_map(0, "")
                cmd.get_wizard().set_map(1, "")
                cmd.get_wizard().set_map(2, "")
                cmd.get_wizard().set_track(1)
                maps = [map1, map2, map3]
                group_name = code.upper() + "_maps_g"
                try:
                    cmd.group(group_name, map1[2])
                    cmd.group(group_name, map2[2])
                except:
                    pass
                for map in maps:
                    cmd.get_wizard().set_map(map[1], map[2])
                    cmd.get_wizard().set_level(map[1], map[3])
                    cmd.color(map[4], map[0] + map[2])
                    try:
                        cmd.group(group_name, map[0] + map[2])
                    except: #no group command available
                        pass
                cmd.get_wizard().update_maps()
            except:
                print(
                    "# Density Wizard could not be loaded with current settings.")
            return
        else:
            typ = typ.strip().split()
        for code in PDBcode:
            for source in typ:
                print("# Fetch %s of %s" % (source, code))
                cmd.fetch(code, type=source)
                if source == "pdb1":
                    cmd.split_states(code)
                    cmd.delete(code)
def unbondatom(selection1="sele", selection2="all", selection="all"):
    '''Unbond atoms.
    USAGE: unbondatom(selection1="sele",selection2="all",selection="sele")
    PYMOL USAGE: bni_unbondatom selection1,selection2,selection
             selection1="sele",
                        "metals", metals in selection (sele)
             selection2="all"
             selection="sele" overall selection1 restriction ("all" for no restiction)
    unbonds selection1 (with restriction 'selection') from selection2
    '''
    selection = _remove_brakets(selection, brakets=0)
    selection1 = _remove_brakets(selection1, brakets=0)
    selection2 = _remove_brakets(selection2, brakets=0)
    metals = ["mg", "zn", "ca", "fe", "mn", "cd", "cu", "co", "ni"]
    info = []
    info.append("# BNI-Tools v. %s" % __version__)
    if selection1 == "metals":
        metal = "+".join(metals)
        selection1 = "symbol %s and (%s)" % (metal, selection)
        try:
            info.append("# Metals:")
            for metal in metals:
                cmd.select("_bni_metal_sele_%s" %
                           metal, "symbol %s and (%s)" % (metal, selection))
                metalcount = cmd.count_atoms("_bni_metal_sele_%s" % metal)
                if metalcount > 0:
                    info.append("#      %s%s:%s metal ions found in selection (%s)." % (
                        metal[0].upper(), metal[1].lower(), metalcount, selection))
            cmd.select("_bni_selection1", selection1)
            atoms1 = cmd.count_atoms("_bni_selection1")
            if atoms1 == 0:
                info.append("# No Metals found in selection (%s)." % selection)
            info.append("# %d Metals selected." % atoms1)
        except:
            info.append("# No atoms selected in (%s)" % selection)
            print("\n".join(info))
            return
    try:
        cmd.select("_bni_selection1", selection1)
        atoms1 = cmd.count_atoms("_bni_selection1")
        if atoms1 == 0:
            info.append("# No atoms in selection. (%s)" % selection)
        info.append("# Unbond %d atoms" % atoms1)
        cmd.unbond(selection1, selection2)
    except:
        info.append("# No atoms selected in (%s)" % selection)
        print("\n".join(info))
        return
    print("\n".join(info))
def combine_sele(selections=None, check="enabled", name="sele"):
    '''Combine selections. 
    USAGE: combine_sele(selections=None,check="enabled",name="sele")
             selections: "sele1 sele2 sele3": add additinal selections
                  check:           "enabled": add only enabled selections
                                       "all": add all  selections
                                  "disabled": add only disabled selectins
                                          -1: do not look for selections
                   name:                sele: specify name of combined selection
    PYMOL: bni_combine_sele selections=None,check=-1,name="sele"
    '''
    selections = _remove_brakets(selections)
    enabled = str(_remove_brakets(check))
    name = str(_remove_brakets(name))
    if enabled == "enabled":# combine enabled selections
        sele_names = set(cmd.get_names("public_selections", enabled_only=1))
    elif enabled == "all":# combine all selections
        sele_names = set(cmd.get_names("public_selections", enabled_only=0))
    elif enabled == "disabled":# combine not enabled selections
        sele1 = set(cmd.get_names("public_selections", enabled_only=1))
        sele2 = set(cmd.get_names("public_selections", enabled_only=0))
        sele_names = sele2 - sele1
    elif enabled == "-1":# do not check only use selections input
        sele_names = []
    sele_names = list(sele_names)
    if selections:
        addsele = selections.split()
        for item in addsele:
            sele_names.append(item)
    if sele_names:
        cmd.select(name, selection="(" + " or ".join(sele_names) + ")")
        #cmd.set_name("_bni"+name,name)
    else:
        cmd.select(name, selection=None)
        #cmd.set_name("_bni"+name,name)
    print("# BNI - Tools v. %s" % __version__)
    print("# Combined selections: %s" % sele_names)
    print("# to (%s)" % name)
def _set_bni_ray(setting="set", app=None, async_=1):
    '''Run sidebar input
    '''
    width = 100
    dpi = 200
    width_unit = "mm"
    if setting == "set" and app:
        value = tkSimpleDialog.askstring('Width mm/cm/inch/dots resolution dpi',
                                         'Please specify width of image and resolution',
                                         parent=app.root)
        value = _remove_brakets(value)
        if value:
            value = value.split()
        else:
            return
        if len(value) == 1:
            width = value[0]
        elif len(value) == 2:# width / dpi
            width = value[0]
            dpi = value[1]
        elif len(value) >= 3:#width width_unit
            width = value[0]
            width_unit = value[1]
            dpi = value[2]
        else:
            print("# Wrong input. Set at least the width.")
            return
        try:
            float(width)
        except ValueError:
            print("# Width has to be a number.")
            return
        units = ["cm", "mm", "inch", "dots"]
        if not width_unit.lower().strip() in units:
            print("# Unit has to be either %s " % (" or ".join(units)))
            return
        try:
            int(dpi)
        except ValueError:
            print("# DPI has to be an integer number.")
            return
    elif setting == "poster":
        print("# POSTER SETTING")
        width = 30
        width_unit = "cm"
        dpi = 300
        print("# width: %s %s" % (width, width_unit))
        print("# resolution:  %s dpi" % dpi)
    elif setting == "journal single":
        print("# JOURNAL SETTING single column")
        width = 86
        width_unit = "mm"
        dpi = 400
        print("# width: %s %s" % (width, width_unit))
        print("# resolution:  %sdpi" % dpi)
    elif setting == "journal double":
        print("# JOURNAL SETTING double column")
        width = 178
        width_unit = "mm"
        dpi = 400
        print("# width: %s %s" % (width, width_unit))
        print("# resolution: %s dpi" % dpi)
    elif setting == "presentation":
        print("# PRESENTATION SETTING")
        width = 150
        width_unit = "mm"
        dpi = 150
        print("# width: %s %s" % (width, width_unit))
        print("# resolution: %s dpi" % dpi)
    elif setting == "hires":
        print("# HIRES SETTING")
        width = 120
        width_unit = "mm"
        dpi = 1200
        print("# width: %s %s" % (width, width_unit))
        print("# resolution: %s dpi" % dpi)
    if async_ == 1:
        print("# Ray is working in the background")
    bni_ray(width, width_unit=width_unit, dpi=dpi, **{'async': async_})
def bni_ray(width, name_png=None, dpi="300", width_unit="cm", **arg):
    '''Ray specified with dpi and width.
    USAGE: bni_ray(width,dpi="300",width_unit="cm",name_png=None,**arg)
            width: width of ray picture in defined units
            name_png:     ./name.png: write png file with ./name.png
                    :               : do not write png just do ray
            width_unit:    : unit of width
                         cm: centimeter  
                       inch: inch
                         mm: millimeter
                       dots: (normal behaviour of ray)
            dpi:        300: dots per inch (resolution)
            **arg: additinal arguments for ray
                   height: integer :{default: 0 (current)}
                antialias: integer :{default: -1 (use antialias setting)}
                    angle: float   :y-axis rotation for stereo image generation
                                    {default: 0.0}
                    shift: float   : x-axis translation for stereo image generation
                                    {default: 0.0}
                 renderer: -1      :respectively
                            0      :default
                            1      :built-in
                            2      :pov-ray, or dry-run {default: 0}
                    async:  0      : do not run in background thread
                            1      : background thread
    PYMOL: bni_ray width,name_png=./my.png,dpi="300",width_unit="cm",**arg e.g. angle=3.0
    '''
    width = float(_remove_brakets(width))
    width_unit = _remove_brakets(width_unit)
    dpi = int(_remove_brakets(dpi))
    if width_unit == "cm":
        calc_res = round(width * dpi / 2.54)
    elif width_unit == "mm":
        calc_res = round(width * dpi / 25.4)
    elif width_unit == "inch":
        calc_res = round(width * dpi)
    elif width_unit == "dots":
        calc_res = int(width)
    else:
        print("Unit not recognized: %s" % width_unit)
        return
    print("# BNI - Tools v. %s" % __version__)
    print("# %d (dots in x direction)" % (calc_res))
    print("# Size(x): %s %s, resolution: %sdpi" % (width, width_unit, dpi))
    print("# [Size(y): calculated automatically.]")
    cmd.ray(width=calc_res, **arg)
    if name_png is not None:
        name_png = _remove_brakets(name_png)
        name = os.path.abspath(os.path.normpath(name_png))
        cmd.png(name, dpi=dpi)
        print("# saved as: %s " % name)
def his_handle(fromto="his2hid", selection="(sele)", hisH_dic=None):
    '''Alter HIS.
    fromto:"his2hid":  HIS to HIP HID HIE
          :"hid2his":  HIS HIP HIE to HIS
    hisH_dic: dictionary of his hydrogen naming convention
              e.g. {"HID":["HD1","H0X"],
                    "HIE":["HE2","H01"]}
    selection: (sele):
    '''
    selection = _remove_brakets(selection, brakets=0)
    fromto = _remove_brakets(fromto)
    try:
        cmd.count_atoms(selection)
    except:
        print("Please select at least one residue in (sele)")
        return
    if fromto == "his2hid":
        print("# BNI - Tools v. %s" % __version__)
        print("# ---------------------")
        print("# Rename HIS to HID HIE HIP")
        print("# Hydrogens have to be added accordingly")
        print("# at the HIS residues ")
        print("# HD1 ?         --> HID")
        print("# HE2 ?         --> HIE")
        print("# HD1 and HE2 ? --> HIP")
        print("# ---------------------")
        allhis = cmd.count_atoms("(r. HIS) and n. CA")
        #hisH_dic naming convention   [reduce,pymol]
        if not hisH_dic:
            hisH_dic = {"HID": ["HD1"],
                        "HIE": ["HE2", "H01"]}
        cmd.select("HIP", "(byres ((n. %s) and r. HIS) and byres ((n. %s) and r. HIS))" % (
            ",".join(hisH_dic.get("HID")), ",".join(hisH_dic.get("HIE"))))
        cmd.alter("(HIP)", "resn='HIP'")
        cmd.sort()
        cmd.select("HID", "(byres ((n. %s) and r. HIS))" %
                   (",".join(hisH_dic.get("HID"))))
        cmd.alter("(HID)", "resn='HID'")
        cmd.sort()
        cmd.select("HIE", "(byres ((n. %s) and r. HIS))" %
                   (",".join(hisH_dic.get("HIE"))))
        cmd.alter("(HIE)", "resn='HIE'")
        cmd.sort()
        print("# %s HIS residues" % allhis)
        print("# -------------------------------")
        print("# %s HID" % cmd.count_atoms("(r. HID and n. CA)"))
        print("# %s HIE" % cmd.count_atoms("(r. HIE and n. CA)"))
        print("# %s HIP" % cmd.count_atoms("(r. HIP and n. CA)"))
        print("# -------------------------------")
        print("# %s HIS residues altered" %
              cmd.count_atoms("((r. HID,HIE,HIP) and n. CA)"))
        print("# -------------------------------")
    elif fromto == "hid2his":
        print("# BNI - Tools v. %s" % __version__)
        print("# --------------------")
        print("# Rename HID HIE HIP to HIS")
        print("# --------------------")
        allhidhiehis = cmd.count_atoms("((r. HID,HIE,HIP) and n. CA)")
        print("# %s HID" % cmd.count_atoms("(r. HID and n. CA)"))
        print("# %s HIE" % cmd.count_atoms("(r. HIE and n. CA)"))
        print("# %s HIP" % cmd.count_atoms("(r. HIP and n. CA)"))
        cmd.select("HIS", "(byres (r. HIE,HIP,HID))")
        cmd.alter("(HIS)", "resn='HIS'")
        cmd.sort()
        allhis = cmd.count_atoms("(r. HIS) and n. CA")
        print("# %s HIS residues" % allhis)
        print("# -------------------------------")
        print("# %s HID,HIE,HIP residues altered" % allhidhiehis)
        cmd.sort()
def polyala(selection="(sele)", name=None):
    '''Get poly alanine chain of selected residues.
    USAGE:  polyala(selection=None,name=None)
              selection: selection of residues
              name:      name of created polyalanine
                         if no name is given the original input
                         will be altered to polyalanine
    PYMOL:  bni_get_polyala (selection),[name]
              selection: selection of residues
              name:      name of created polyalanine
                         if no name is given the original input
                         will be altered to polyalanine
    ADDITIONAL: if "name" is specified the structure looks wired because of the
                lack of bonds between the altered alanine residues.
                If you want to see the bonds you can save and reopen the struture.
                (Pymol recognize the new bonds after loading)
    '''
    if selection is None:
        print("# Please select at least one residue.")
        return
    tmpname2 = "_polyala2"
    tmpname3 = "_polyala3"
    selection = _remove_brakets(selection)
    try:
        cmd.count_atoms(selection)
    except:
        print("# Empty selection. Select at least one atom.")
        return
    if name != None:
        if name in cmd.get_names():
            cmd.delete(name)
        cmd.select(tmpname2, "((byobj %s) and not (byres %s))" %
                   (selection, selection))
        cmd.create(
            tmpname3, "((byres %s) and (n. CB,O,CA,N,C) and not hetatm and not solvent)" % (selection))
        cmd.alter(("%s and not r. GLY") % tmpname3, "resn='ALA'")
        cmd.sort()
        cmd.create(name, ("(%s or %s)" % (tmpname3, tmpname2)))
        cmd.delete(tmpname2)
        cmd.delete(tmpname3)
    else:
        cmd.remove(
            "((byres %s) and not (n. CB,O,CA,N,C) and not hetatm and not solvent)" % selection)
        cmd.select(
            tmpname3, "((byres %s) and (n. CB,O,CA,N,C) and not hetatm and not solvent)" % (selection))
        cmd.alter(("%s and not r. GLY") % tmpname3, "resn='ALA'")
        cmd.sort()
        cmd.delete(tmpname3)
def cgo_point(p):
    cgo = cgo_dic()
    COLOR = cgo["COLOR"]
    SPHERE = cgo["SPHERE"]
    x, y, z = p
    return [SPHERE, float(x), float(y), float(z), 0.2]
def get_triangle(selection, name="tri", invert=1, sphere=0):
    '''Get a triangle from 3 atom selection.
    USAGE: get_triangle(selection,name="tri",size=100,invert=0,sphere=1)
             name     :        :name of the triangle
             selection:        :selection containing three atoms
             invert   : int    : 2: invert normal for lightning
                                 1: do not invert normal
                                 0: do not create normal
             sphere   : int    : 1: draw spheres on edges of the plane
                               : 0: do not draw spheres
    PYMOL: bni_get_triangle selection,name,invert,sphere
    '''
    name = _remove_brakets(name)
    selection = _remove_brakets(selection, brakets=0)
    invert = _remove_brakets(invert)
    invert = int(invert)
    coordinates = []
    try:
        mol = cmd.get_model(selection)
        mol.sort()
    except:
        print("# Empty selection %s" % selection)
        return
    else:
        if len(mol.atom) != 3:
            print("# Please select exactly three atoms.")
            return
        for a in mol.atom:
            coordinates.append(a.coord)
        tri = get_triangle_coord(
            points=coordinates, invert=invert, sphere=sphere)
        triName = check_names(name)
        cmd.load_cgo(tri, triName)
        cmd.show("cgo", triName)
def get_triangle_coord(points=((0, 0, 0), (0, 1, 0), (1, 0, 0)), invert=0, sphere=1):
    '''Get coords for triangle polygon.
    Currently only 3 points are supported.
    '''
    #invert=0
    cgo = cgo_dic()
    END = cgo["END"]
    obj = []
    rounds = len(points) - 1
    centre = getCenterCoord(points)
    # create triangles from midpoint of points (useful for later versions with
    # more points)
    point_index = [0, 1, 2]
    for i in range(rounds):
        #cpv.normalize(cpv.sub(points[i], centre))
        vector1 = cpv.normalize(cpv.sub(centre, points[i]))
        #cpv.normalize(cpv.sub(points[i+1], centre))
        vector2 = cpv.normalize(cpv.sub(centre, points[i + 1]))
        #vector3 = cpv.normalize(cpv.sub(points[i], points[i+1]))
        #if normal_indication==1:
        #normal = cpv.normalize(cpv.cross_product(vector1, vector2))
        normal = cpv.cross_product(vector2, vector1)
        #
        obj_new = triangle([centre, points[i], points[i + 1]], norm=normal,
                           invert=invert, sphere=sphere, point_index=point_index, finishobj=0)
        obj.extend(obj_new)
    if rounds > 1:
        vector1 = cpv.normalize(cpv.sub(centre, points[-1]))
        vector2 = cpv.normalize(cpv.sub(centre, points[0]))
        normal = cpv.cross_product(vector2, vector1)
        obj_new = triangle([centre, points[-1], points[0]], norm=normal,
                           invert=invert, sphere=sphere, point_index=point_index, finishobj=0)
        obj.extend(obj_new)
    obj.append(END)
    return obj
def triangle(points=[(0, 0, 0), (1, 1, 0), (1, 0, 0), (0, 1, 0)], norm=None, invert=0, sphere=1, point_index=[0, 1, 2, 0, 2, 3], finishobj=1):
    '''Get triangles.
          points: unique list of points for triangles in same plane.
          norm: normal of triangles in plane
          invert: 0,1 : invert normal
          sphere: 0,1 :show additional spheres on edges
          point_index: list of point indices to create triangles
          finishobj: 0,1 : add END to object or not
    '''
    cgo = cgo_dic()
    BEGIN = cgo["BEGIN"]
    COLOR = cgo["COLOR"]
    TRIANGLE_STRIP = cgo["TRIANGLE_STRIP"]
    TRIANGLES = cgo["TRIANGLES"]
    NORMAL = cgo["NORMAL"]
    VERTEX = cgo["VERTEX"]
    END = cgo["END"]
    obj = []
    if sphere == 1:
        for p in points:
            obj.extend(cgo_point(p))
    obj.extend([BEGIN, TRIANGLES])
    corners = []
    for indexes in point_index:
        corners.append(points[indexes])
    if invert == 2 and norm is not None:
        normal = norm
    elif invert == 1 and norm is not None:
        normal = [i * (-1) for i in norm]
    else:  
        for corner in corners:
            obj.append(VERTEX)
            obj.extend(corner)
        if finishobj == 1:
            obj.append(END)
        return obj
    for corner in corners:
        obj.append(NORMAL)
        obj.extend(normal)
        obj.append(VERTEX)
        obj.extend(corner)
    if finishobj == 1:
        obj.append(END)
    return obj
def get_plane_coord(p1, p2, p3, cushion, invert=0, sphere=1):
    '''Get coords of plane.
    '''
    vector1 = cpv.normalize(cpv.sub(p2, p1))
    vector2 = cpv.normalize(cpv.sub(p3, p1))
    normal = cpv.cross_product(vector1, vector2)
    centre = getCenterCoord([p1, p2, p3])
    vector2 = cpv.cross_product(normal, vector1)
    x = cpv.scale(vector1, cushion)
    y = cpv.scale(vector2, cushion)
    c1 = cpv.add(cpv.add(centre, x), y)
    c2 = cpv.sub(cpv.add(centre, x), y)
    c3 = cpv.sub(cpv.sub(centre, x), y)
    c4 = cpv.add(cpv.sub(centre, x), y)
    return triangle([c1, c2, c3, c4], norm=normal, invert=invert, sphere=sphere, point_index=[0, 1, 2, 0, 2, 3])
def get_plane_handler(selection="(sele)", app=None, action="set", name="plane", size=100, invert=2, sphere=0):
    '''Get a plane using get_plane GUI.
    '''
    action = _remove_brakets(action)
    name = check_names(name)
    size = size
    invert = str(invert)
    sphere = str(sphere)
    if action == "set" and app:
        value = tkSimpleDialog.askstring('Plane settings',
                                         'Specify [size %s] [invert 2/1/0]  [sphere 1/0]' % "A",
                                         parent=app.root)
        value = _remove_brakets(value)
        if value:
            value = value.split()
        else:
            return
        if len(value) == 1:
            size = value[0]
        elif len(value) == 2:
            size = value[0]
            invert = str(value[1])
        elif len(value) >= 3:
            size = value[0]
            invert = str(value[1])
            sphere = str(value[2])
        else:
            print("# Wrong input. Set at least the size.")
            return
        try:
            float(size)
        except ValueError:
            print("# Size has to be a number.")
            return
        units = ["0", "1", "2"]
        if not invert.lower().strip() in units:
            print("# Invert has to be either %s " % (" or ".join(units)))
            return
        units = ["1", "0"]
        if not sphere.lower().strip() in units:
            print("# Sphere has to be either %s " % (" or ".join(units)))
            return
    elif action == "small":
        size = 5
        invert = 1
        sphere = 0
    elif action == "medium":
        size = 10
        invert = 1
        sphere = 0
    elif action == "large":
        size = 20
        invert = 1
        sphere = 0
    get_plane(selection, name=name, size=size, invert=invert, sphere=sphere)
def get_plane(selection, name="plane", size=100, invert=0, sphere=1):
    '''Get a plane from 3 atom selection.
    USAGE: get_plane(selection,name="plane",size=100,invert=0,sphere=1)
             name     :        :name of the plane
             selection:        :selection containing three atoms
             size     : int    :cushion size of the plane around selection
             invert   : int    : 2: invert normal for lightning
                                 1: do not invert normal
                                 0: do not create normal
             sphere   : int    : 1: draw spheres on edges of the plane
                               : 0: do not draw spheres
    PYMOL: bni_get_plane selection,name,size,invert,sphere
    '''
    name = _remove_brakets(name)
    selection = _remove_brakets(selection, brakets=0)
    size = _remove_brakets(size)
    invert = _remove_brakets(invert)
    size = int(size)
    invert = int(invert)
    coordinates = []
    try:
        mol = cmd.get_model(selection)
    except:
        print("# Empty selection %s" % selection)
        return
    else:
        if len(mol.atom) != 3:
            print("# Please select exactly three atoms.")
            return
        for a in mol.atom:
            coordinates.append(a.coord)
        p1 = coordinates[0]
        p2 = coordinates[1]
        p3 = coordinates[2]
        plane = get_plane_coord(p1, p2, p3, int(size), invert, sphere)
        planeName = name
        cmd.load_cgo(plane, planeName)
        cmd.show("cgo", planeName)
def _get_list(input, message="Error"):
    '''Remove brakets and return list.
       Returns a list from input,
       else a error massage is printed
    '''
    input = _remove_brakets(input)
    if type(input) is str:
        try:
            input = map(int, input.split(","))
        except TypeError:
            print("# %s:\n%s" % (message, mobile))
            return
        except ValueError: 
            input = input.split(",")
    return list(input)
def _maxnamelenght(types):
    '''Get max lenght of names.
    types: 1 : get name lenght
           2 : get name and title lenght
    '''
    stdwidth = sidebar_stdwidth
    multipl = 12.0 
    addpixel = 1 * multipl
    if __pymol_version__ > 1.4:
        addpixel = 6 * multipl
    names = cmd.get_names("public_selections")
    names.extend(cmd.get_names("public_objects"))
    if types == 1:
        try:
            maxlenght = max([len(i) for i in names])
        except ValueError:
            return stdwidth
        try:
            namelenght = int(
                round(float(maxlenght * multipl + addpixel + 0.5)))
        except ValueError:
            return stdwidth
        if namelenght < stdwidth:
            namelenght = stdwidth
        return namelenght
    elif types == 2:
        titles = []
        for tit in names:
            try:
                title = cmd.get_title(tit, cmd.get_state())
                if title == None:
                    title = ""
            except AttributeError:
                titles.append("")
            titles.append(title)
        if len(names) != len(titles):
            print("# WARNING: name and title mismatch/nnames:%s/ntitles:%s" %
                  (names, titles))
        else:
            try:
                maxlenght = max([len(u) + len(v)
                                 for u, v in zip(names, titles)])
                #maxlenght=maxlenght+2*multipl
            except ValueError:
                return stdwidth
        try:
            namelenght = int(
                round(float(maxlenght * multipl + addpixel + 0.5)))
        except ValueError:
            return stdwidth
        if namelenght < stdwidth:
            namelenght = stdwidth
        return namelenght
def runsidebar(set="", value=None, app=None):
    '''Run sidebar input
    '''
    if value == "set" and app:
        value = tkSimpleDialog.askstring('Input Sidebar Value',
                                         'Please specify a Value e.g 250',
                                         parent=app.root)
    sidebar(set="width", value=value)
def sidebar(set="", value=None):
    '''Sidebar Settings.
    USAGE: sidebar(set,value)
       set: 
              "on"   :  set sidebar on
              "off"  :  set sidebar off
              "width":  set sidebar to value
              "size" :  set hight of each sidebar entry
              "overlayon": sidebar overlay on
              "overlayoff": sidebar overlay off
       value:
              value : int : value of sidebar
    '''
    valuedic = {}
    if set == "on":
        cmd.set("internal_gui", "on")
    elif set == "off":
        cmd.set("internal_gui", "off")
    elif set == "overlayon":
        try:
            if __pymol_version__ < 1.2:
                cmd.set("internal_gui_mode", 1)
            else:
                cmd.set("internal_gui_mode", 2)
        except:
            pass
    elif set == "overlayoff":
        try:
            cmd.set("internal_gui_mode", 0)
        except:
            pass
    elif set == "width":
        if type(value) != type(None):
            if value == "fit1":#fit to names
                fit = _maxnamelenght(1)
                valuedic.update({"width": fit})
            elif value == "fit2":#fit to names and titles
                fit = _maxnamelenght(2)
                valuedic.update({"width": fit})
            else:
                try:
                    value = int(float(value))
                    valuedic.update({"width": value})
                except:
                    print("# Wrong value:", value)
                    return
        else:
            return
        try:
            cmd.set("internal_gui_width", valuedic.get("width", 250))
            print("# Set sidebar to %s" % valuedic.get("width", 250))
        except:
            return
    elif set == "size":
        cmd.set("internal_gui_control_size", valuedic.get("size", 18))
def _check_selection(selection="(sele)"):
    '''Check if selection contains atoms.
    '''
    try:
        return cmd.count_atoms(selection)
    except:
        return 0
def get_surface_handler(selection="(sele)", action="mse", app=None):
    '''Surface representation GUI using mse.
    '''
    selection = _remove_brakets(selection)
    action = _remove_brakets(action)
    len_atoms = _check_selection(selection)
    if action == "Standard":
        mse(sele=selection + " and (hetatm or solvent)",
            stype="-1", ssurface="0")
        mse(sele=selection + " and ! (hetatm or solvent)",
            stype="-1", ssurface="1")
        print("# BIN - Tools v. %s" % __version__)
        print("# surface representation flags are set to standard values")
        print("# for ( %s ). (current present objects only)." % (selection))
        print("# show the surface to see the effect")
        return
    if action == "Reset":
        mse(sele="(hetatm or solvent)", stype="-1", ssurface="0")
        mse(sele="!(hetatm and solvent)", stype="-1", ssurface="1")
        print("# BIN - Tools v. %s" % __version__)
        print("# all surface representation flags are set to standard values")
        print("# (current present objects only)." % (selection))
        print("# show the surface to see the effect")
        return
    if action == "Methionine" or action == "mse":
        mse(sele=selection + " and r. MSE", stype="-1", ssurface="1")
        print("# BIN - Tools v. %s" % __version__)
        print("# Methionine surface representation is set to on")
        print("# for ( %s ). (current present objects only)." % (selection))
        print("# show the surface to see the effect")
        return
    if len_atoms == 0:
        print("# No atoms in selection %s. Select at least one atom." %
              selection)
        return
    elif action == "Include":
        mse(sele=selection, stype="-1", ssurface="1")
        print("# BIN - Tools v. %s" % __version__)
        print("# Included %s atom(s) to surface calculation." % len_atoms)
        print("# show the surface to see the effect")
    elif action == "Exclude":
        mse(sele=selection, stype="-1", ssurface="0")
        print("# BIN - Tools v. %s" % __version__)
        print("# Excluded %s atom(s) from surface calculation." % len_atoms)
        print("# show the surface to see the effect")
    elif action == "Show":
        mse(sele=selection, stype="-1", ssurface="1")
        print("# BIN - Tools v. %s" % __version__)
        print("# Show surface for %s atom(s)." % len_atoms)
        print("# show the surface to see the effect")
    elif action == "Hide":
        mse(sele=selection, stype="-1", ssurface="2")
        print("# BIN - Tools v. %s" % __version__)
        print("# Hide surface for %s atom(s)." % len_atoms)
        print("# show the surface to see the effect")
def delalternate(sele="sele"):
    '''Delete alternates.
    '''
    sele = _remove_brakets(sele, brakets=0)
    selection = "(" + sele + " and not (alt '' or alt A))"
    try:
        len_atoms = cmd.count_atoms(selection)
        len_residues = cmd.count_atoms(
            "(byres " + selection + ") and n. CA and not e. CA and (alt '' or alt A)")
    except:
        print("# No residues in selection to delete.")
        return
    cmd.remove(selection)
    cmd.sort()
    print("# BIN - Tools v. %s" % __version__)
    print("# Deleted %d alternate atoms in %d residues." %
          (len_atoms, len_residues))
def altermse(sele="(sele)"):
    '''Edit MSE to MET permanantly.
    '''
    sele = _remove_brakets(sele, brakets=0)
    try:
        cmd.count_atoms(sele)
    except:
        print("# No atoms in selection %s" % sele)
        return
    selection = "(" + sele + " and r. MSE and n. SE)"
    cmd.alter(selection=selection, expression="name='SD'", quiet=1, space=None)
    cmd.sort()
    selection = "(" + sele + " and r. MSE and n. SD)"
    cmd.alter(selection=selection, expression="elem='S'", quiet=1, space=None)
    cmd.sort()
    cmd.flag(25, "(" + sele + " and r. MSE)", "clear")
    cmd.sort()
    selection = "(" + sele + " and r. MSE)"
    mses = cmd.count_atoms(selection + " and n. CA and not e. CA")
    cmd.alter(
        selection=selection, expression="type='ATOM'", quiet=1, space=None)
    cmd.sort()
    cmd.alter(
        selection=selection, expression="resn='MET'", quiet=1, space=None)
    print("# BIN - Tools v. %s" % __version__)
    print("# Changed %d MSE residues to MET." % mses)
    cmd.sort()
def mse(sele="r. MSE", stype="-1", ssurface="1"):
    '''Handle MSE atom surface representation (Seleno- Methionine).
    USAGE:  bni_mse sele="r. MSE" stype=-1 ssurface=1
            sele:     selection : standard is "r. MSE"
            stype:    1: set MSE type from HETATM to ATOM
                      0: set MSE type from ATOM to HETATM
                     -1: do nothing (standard)
            ssurface: 1: set surface representation to 1 for MSE (create surface) (standard)
                      2: set surface representation to 2 for MSE (dont show but include)
                      0: set surface representation to 0 for MSE (exclude from surface)
                         this is standard within pymol
                     -1: do nothing
    '''
    sele = _remove_brakets(sele, brakets=0)
    if stype == "1":
        cmd.alter("%s" % sele, 'type="ATOM"')
    elif stype == "0":
        cmd.alter("%s" % sele, 'type="HETATM"')
    else:
        pass
    if ssurface == "1":
        cmd.flag(25, sele, "clear")
        cmd.flag(24, sele, "clear")
    elif ssurface == "2":
        cmd.flag(25, sele, "clear")
        cmd.flag(24, sele, "set")
    elif ssurface == "0":
        cmd.flag(25, sele, "set")
        cmd.flag(24, sele, "clear")
    else:
        pass
    cmd.sort()
def _redirect_alter(selection, property, addprint=1, command="alter"):
    '''redirect standart out temporarily to get cmd.alter command variables.
     selection:
     property: index
     addprint:1: additional print to stdout
     command:"alter": use the cmd.alter command
             "iterate": use the cmd.iterate command (v.1.2r1 also needs print redirection!)
    '''
    if command == "alter":
        command = cmd.alter
    elif command == "iterate":
        command = cmd.iterate
    values = []
    class MyRedirect:
        def __init__(self, stdout):
            self.stdout = stdout
        def write(self, s):
            if s.strip():
                values.append(s)
            if addprint == 1:
                self.stdout.write(s)
    old_stdout = sys.stdout
    sys.stdout = MyRedirect(sys.stdout)
    # index of the atoms is now in a list(values) as string
    command(selection, property)
    sys.stdout = old_stdout
    return values
def _seq_format(selection="sele", wrap=60, format="fasta", gap=1, gap_sign="-"):
    '''sequence output.
    format: pir
            modeller
            fasta
    wrap: sequence chars in one line
    gap:      1:
              0:
    gap_sign: -:
              .:
              _:
          lower:gap AA are written in lower case
    returns (fromated sequences,gap_residues)
    '''
    output_format = "'%s;%s'%(resn,resi)"
    chains = cmd.get_chains(selection)
    allnu = []
    try:
        names = cmd.get_names("public_objects")
    except:
        names = cmd.get_names()
    exc = HP()
    nr_out = []
    for obj in names:
        for chain in chains:
            if chain == '':
                chain = "''"
            select_new = "%s and " % obj + selection + " and c. %s" % chain
            u = _redirect_alter(
                select_new, output_format, addprint=0, command="alter")
            #print select_new
            if not u:
                continue
            u = [((_remove_brakets(i.split(";")[0]), _remove_brakets(i.split(";")[1])))
                 for i in u if type(i) is type("string")]
            lenght = len(u)
            aa = [a[0] for a in u]
            nr = [a[1] for a in u]
            aa = exc.three2one(aa)
            if not len(aa) == len(nr):
                print("# WARNING: AA and NR do not have the same lenght.")
                print(
                    "# WARNING: output in line 2 of pir format may be incorrect.")
            nu = "".join(aa)
            if chain == '':
                chain = "."
            #compare with whole chain to insert gaps as "-"
            if gap == 1:
                select_gap = "(bychain " + select_new + \
                    ") and n. CA and not e. CA and (alt '' or alt A)"
                u_gap = _redirect_alter(select_gap, output_format, addprint=0)
                u_gap = [((_remove_brakets(i.split(";")[0]), _remove_brakets(
                    i.split(";")[1]))) for i in u_gap if type(i) is type("string")]
                lenght_gap = len(u_gap)
                aa_gap = [a_gap[0] for a_gap in u_gap]
                nr_gap = [a_gap[1] for a_gap in u_gap]
                aa_gap = exc.three2one(aa_gap)
                if not len(aa_gap) == len(nr_gap):
                    print(
                        "# WARNING: AA_gap and NR_gap do not have the same lenght.")
                    print("# WARNING: output may be incorrect.")
                #nu_gap="".join(aa_gap)
                #compare and create new seq with gaps
                new_seq = []
                #print nr_gap,aa_gap
                #print nr,aa
                nr_out.append("# Residue from - to")
                nr_out.append("# for %s:%s:%d (out of %d)" %
                              (obj, chain, lenght, len(nr_gap)))
                contains_gap = 0
                new_nr = []
                for number, amino in zip(nr_gap, aa_gap):
                    if number in nr:
                        new_seq.append(amino)
                        new_nr.append((number, amino))
                    else:
                        if new_nr:
                            nr_out.append("%1s %4s - %1s %4s : (%3d AA)" % (new_nr[0][1], new_nr[0][
                                          0], new_nr[-1][1], new_nr[-1][0], len(new_nr))) 
                            contains_gap = 1
                            new_nr = []
                        if gap_sign == "lower":
                            new_seq.append(amino.lower())
                        else:
                            new_seq.append(gap_sign)
                if contains_gap == 0:# no gaps
                    nr_out.append("%1s %4s - %1s %4s : (%3d AA)" %
                                  (aa_gap[0], nr_gap[0], aa_gap[-1], nr_gap[-1], len(nr_gap)))
                elif new_nr: 
                    nr_out.append("%1s %4s - %1s %4s : (%3d AA)" % (new_nr[0][1], new_nr[0][
                                  0], new_nr[-1][1], new_nr[-1][0], len(new_nr))) 
                    new_nr = []
                    if gap_sign == "lower":
                        new_seq.append(amino.lower())
                    else:
                        new_seq.append(gap_sign)
                nu = "".join(new_seq)
            if format == "modeller":
                allnu.append(">P1;%s" % (obj))
                label = "structureX"
                source = "source"
                rfac = -1.0
                resolution = -1.0
                if gap != 1:
                    allnu.append("%s:%s:%s:%d:%s:%d:SeqLenght %d:%s:%.2f:%.2f" % (
                        label, obj, chain, int(nr[0]), chain, int(nr[-1]), lenght, source, rfac, resolution))
                else:
                    allnu.append("%s:%s:%s:%d:%s:%d:SeqLenght %d:%s:%.2f:%.2f" % (label, obj, chain, int(
                        nr_gap[0]), chain, int(nr_gap[-1]), lenght, source, rfac, resolution))
                #every "wrap" res return
                nu = _seq_wrap(nu, wrap)
                nu[-1] += "*"
            elif format == "pir":
                allnu.append(">P1;%s" % (obj))
                allnu.append("")
                #every "wrap" res return
                nu = _seq_wrap(nu, wrap)
                nu[-1] += "*"
            elif format == "fasta":
                if gap == 1:
                    allnu.append(">%s:%s:%d out of %d:" %
                                 (obj, chain, lenght, len(nr_gap)))
                else:
                    allnu.append(">%s:%s:%d" % (obj, chain, lenght))
                #every "wrap" res return
                nu = _seq_wrap(nu, wrap)
            allnu.extend(nu)
    return allnu, nr_out
def _seq_wrap(seq, wrap=60):
    '''Sequence wrapper.
    Wraps a sequence into a list.
    seq: string input
    wrap: 60: chars per list entry
    '''
    slice_begin = 0
    slice_end = wrap
    seq_list = []
    while 1:
        seq_slice = seq[slice_begin:slice_end]
        if len(seq_slice) < wrap:
            seq_list.append(seq_slice)
            break
        seq_list.append(seq_slice)
        slice_begin = slice_end
        slice_end += wrap
    return seq_list
def _mw(values=[]):
    '''Return MW  and variance and len.
    '''
    sum1 = sum(values)
    num = len(values)
    try:
        mean = float(sum1) / float(num)
        std = (sum([(a - mean) ** 2 for a in values]) / (num)) ** 0.5
    except ZeroDivisionError:
        mean = 0.
        std = 0.
    return mean, std, num
def get_sequence(selection="sele", selecting="residue", format="normal", name=None, wrap="60", gap="1", gap_sign="-"):
    '''Print selected residues.
         USAGE: get_sequence(selection="sele",selecting="residue",format="normal",name=None,wrap="60")
                selection: sele    : selection of residues or atoms.
                selecting: residue: get all residues of selection (without alternates)
                              atom: get all atoms of selection
                             chain: get all chains of selection
                    name: None    : name of selection for output
                   format: normal: object: chain: residue-residue number
                          normalb: calc bfac per residue
                     normalb_main: calc bfac per residue (only mainchain)
                     normalb_side: calc bfac per residue (only sidechain)
                          normalq: calc occ per residue
                     normalq_main: calc occ per residue (only mainchain)
                     normalq_side: calc occ per residue (only mainchain)
                            fasta: amino acid one letter code in fasta format
                              pir:    amino acid one letter code in pir format
                         modeller:    amino acid one letter code in Modeller pir format
                   wrap:  60: for sequence output wrap at 60 residues
                   gap:  1: add gaps between selected amino acids
                         0: do not add gaps
                   gap_sign: -   : char of gap 
                                  "lower": gap amino acids are written as lower case
         PYMOL:  bni_get_sequence selection=sele,selecting=residue,format=normal,wrap=60,gap_sign=lower
         ADDITIONAL: alternates are not considered! eighter alternate A or alterneta "" are selected
    '''
    format = _remove_brakets(format)
    selecting = _remove_brakets(selecting)
    selection = _remove_brakets(selection)
    wrap = int(_remove_brakets(wrap))
    gap = int(_remove_brakets(gap))
    gap_sign = _remove_brakets(gap_sign)
    try:
        a = cmd.count_atoms(selection)
    except:
        a = None
    if not a:
        print("# Select at least one atom.")
        return
    if name is None:
        name = selection
    prelist = []
    prelist.append("# BNI-Tools v. %s" % __version__)
    prelist.append("# ---------------")
    prelist.append("# Output for %s" % name)
    prelist.append("# ---------------")
    if selection == "enabled":
        pass
    elif selection == "enabled1":
        pass
    else:
        output_format = '"%s:%s: %s-%s"%(model,chain,resn,resi)'
        if selecting == "residue":
            selection = "(byres " + selection + ")"
            selection = "(" + selection + " and (alt '' or alt A))"
            if format == "normal" or format[:7] == "normalb" or format[:7] == "normalq":
                output_format = '"%s:%s: %s-%s "%(model,chain,resn,resi)'
                u = _redirect_alter(selection, output_format, addprint=0)
                u = [_remove_brakets(i) for i in u]
                resis = []
                atoms = []
                for resi in u:
                    r = resi.split(":")
                    #obj:chain:resn:resi
                    radd = (r[0], r[1], r[2].split("-")[0], r[2].split("-")[1])
                    if radd not in resis:
                        atoms.append(resi)
                        resis.append(radd)
                u = atoms[:]
                if format[:7] == "normalb" or format[:7] == "normalq":
                    if format[:7] == "normalb":
                        val = "b"
                        output_format = '"%4.2f"%(b)'
                    else:
                        val = "q"
                        output_format = '"%4.2f"%(q)'
                    resi = []
                    val_allatoms=[]
                    for residue in u:
                        residue = residue.split(":")
                        selection = "(" + residue[0] + " and c. " + residue[1] + " and r. " + residue[
                            2].split("-")[0] + " and i. " + residue[2].split("-")[1] + ")"
                        selection = "(" + selection + \
                            " and (alt '' or alt A) and !n. H)"
                        if format[7:] == "_main":
                            selection = "(" + selection + " and (n. CA,C,N,O) or (" + \
                                selection + ") and !polymer and !n. H)"
                        elif format[7:] == "_side":
                            selection = "(" + selection + " and (!(n. CA,C,O or (n. N and !r. PRO))) or (" + \
                                selection + ") and !polymer and !n. H)"
                        a = _redirect_alter(
                            selection, output_format, addprint=0)
                        a = [_remove_brakets(i) for i in a]
                        val_list = [float(i) for i in a]
                        val0, val1, val2 = _mw(val_list)
                        residue.extend(
                            ["%s= %5.2f" % (val, val0), "%s= %5.2f" % ("s", val1), "atoms %2d" % val2])
                        val_allatoms.extend(val_list)
                        resi.append(":".join(residue))
                    # calc mean for all atoms
                    valall0,valall1,valall2 = _mw(val_allatoms)
                    u = resi
                if format[:7] == "normalb" or format[:7] == "normalq":
                    prelist.append("# No hydrogens and no alternates are considered in calculation")
                if format[7:] == "_main":
                    prelist.append(
                        "# Only mainchain atoms are considered in calculation.")
                elif format[7:] == "_side":
                    prelist.append(
                        "# Only sidechain atoms are considered in calculation.")
                prelist.extend(u)
                if format[:7] == "normalb" or format[:7] == "normalq":
                    prelist.append("#----------")
                    prelist.append("All atoms:%s= %5.2f:%s= %5.2f:atoms %2d" % (val, valall0,"s",valall1,valall2))
                print("\n".join(prelist))
                return u, []
            elif format == "fasta" or format == "pir" or format == "modeller":
                selection = "(" + selection + " and n. CA and not e. CA )"
                (allnu, nr_gap) = _seq_format(selection=selection,
                                              wrap=wrap,
                                              gap=gap,
                                              gap_sign=gap_sign,
                                              format=format)
                if gap == 1:
                    prelist.extend(nr_gap)
                prelist.append("# SEQUENCE FORMAT: %s" % format)
                prelist.append("# ---------------------------")
                prelist.extend(allnu)
                print("\n".join(prelist))
                return allnu, nr_gap
            else:
                output_format = '"%s"%(resn)'
                u = _redirect_alter(selection, output_format, addprint=0)
                #print u
                u = [_remove_brakets(i) for i in u]
                print("\n".join(u))
                return u, []
        elif selecting == "atom":
            selection = selection
            output_format = '"%s:%s: %s-%s :%s"%(model,chain,resn,resi,name)'
            u = _redirect_alter(selection, output_format, addprint=0)
            #print u
            u = [_remove_brakets(i) for i in u]
            print("\n".join(u))
            return u, []
        elif selecting == "chain":
            u = cmd.get_chains(selection)
            print("\n".join(u))
            return u, []
        #cmd.alter(selection,output_format)
class Sym:
    def __init__(self):
        self.org_cc = []  
        self.sym_cc = []  
        self.org_sel = [] 
        self.org_sel_show = []
        self.sym_sel = [] 
        self._cutoff = 0.0
        self.name_sym_obj = ""
        self.picked_selection = ""
        self.ngroup = None
    def create_group(self, gname="bni_sym"):
        '''Create group.
        '''
        self.ngroup = bni_group(gname)
    def groupit(self, names=[]):
        '''Try to group names.
        '''
        if self.ngroup:
            try:
                cmd.group(self.ngroup, " ".join(names))
            except:
                pass
    def _setcutoff(self, restrict=6.0):
        '''Check if cutoff has to be restricted.
        '''
        if self._cutoff > restrict:
            self._cutoff = restrict
    def create_sym(self, cutoff=4.0, selection="sele"):
        '''Create smmetry mates of selection.
        USAGE: 
                create_sym([cutoff=4.0[,selection="sele"]])
        RETURNS: 
                selection_name(s) of symetry involved atoms of symetry mates
        '''
        selection = _remove_brakets(selection)
        self.picked_selection = selection
        self._cutoff = cutoff
        self._setcutoff()
        try:
            self.name_sym_obj = cmd.get_object_list(selection)[0]
        except:
            self.name_sym_obj = ""
        cmd.select("_sym_sele", "(byobj %s)" % selection)
        cmd.create("_sym", "_sym_sele")
        name_symetry_mate = "s_%s_" % self.name_sym_obj
        name_symetry_mate = check_names(name_symetry_mate)
        cmd.symexp(name_symetry_mate, "_sym", selection, cutoff=cutoff, segi=1)
        cmd.delete("_sym_sele")
        cmd.delete("_sym")
        select_name = "_sym_%s_%1.1f" % (self.name_sym_obj, self._cutoff)
        select_name = check_names(select_name)
        #remove
        name_symetry_mates = name_symetry_mate + "*"
        cmd.select(select_name, "((byobj %s) around %1.1f) and not ((byobj %s) and not %s ) and %s" % (
            selection, self._cutoff, selection, selection, name_symetry_mates), enable=0)
        cmd.select(select_name[1:], "(%s around %1.1f) and not ((byobj %s) and not %s ) and %s" % (
            selection, self._cutoff, selection, selection, name_symetry_mates), enable=0)
        symetrys = [r for r in cmd.get_names() if name_symetry_mate in r]
        colors = 2
        self.create_group()
        self.groupit(symetrys)
        self.groupit([select_name[1:]])
        for item in symetrys:
            try:
                util.cba(colors, item)
            except:
                pass
            colors += 1
        return select_name
    def create_org(self, cutoff=4.0, sele=None):
        '''Create selection of orig atoms involved within symmetry.
        USAGE:
                create_org([cutoff=4.0
                           [,sele="s_*"
                           ]])
        '''
        if sele is None:
            sele = "s_%s*" % self.name_sym_obj
        self._cutoff = cutoff
        #restrict selection cutoff to a reasonable distance (max. 6.0 A)
        self._setcutoff()
        cmd.select("_sym", sele)
        #restrict selection cutoff to a reasonable distance (max. 6.0 A)
        select_name = "_org_%s_%1.1f" % (self.name_sym_obj, self._cutoff)
        cmd.select(select_name[1:], "(_sym around %1.1f) and (%s)" % (
            self._cutoff, self.picked_selection), enable=0)
        cmd.select(select_name, "(_sym around %1.1f)" %
                   (self._cutoff), enable=0)
        cmd.delete("_sym")
        self.groupit([select_name[1:]])
        return select_name
    def show_org(self, select_name, show_name="surface", cutoff=4.0, surf_transparency=0.0):
        '''Create different views of symmetry involved selection.
        USAGE: show_org(select_name,
                        [show_name="surface"
                        [,cutoff=4.0,
                        [,surface_transparency=0.0,
                        ]])
                        show_name:"sticks"
                                  "surface"
        '''
        select_name = _remove_brakets(select_name)
        self._cutoff = cutoff
        #restrict selection cutoff to a reasonable distance (max. 6.0 A)
        self._setcutoff()
        cmd.enable(select_name)
        cmd.show(show_name, "_org_%s_%1.1f" %
                 (self.name_sym_obj, self._cutoff))
        if show_name == "surface":
            cmd.set("transparency", surf_transparency)
        cmd.disable("_org_%s_%1.1f" % (self.name_sym_obj, self._cutoff))
    def show_sym(self, select_name, show_name="surface", cutoff=4.0, surf_transparency=0.0, object=None):
        '''Create different views of symmetry involved selection.
        USAGE: show_sym(select_name,
                        [show_name="surface"
                        [,cutoff=4.0,
                        [,surface_transparency=0.0,
                        [,sele="s_*"]]])
                        show_name:"sticks"
                                  "surface"
        '''
        select_name = _remove_brakets(select_name)
        self._cutoff = cutoff
        #restrict selection cutoff to a reasonable distance (max. 6.0 A)
        self._setcutoff()
        cmd.enable(select_name)
        cmd.show(show_name, "_sym_%s_%1.1f" %
                 (self.name_sym_obj, self._cutoff))
        if object is None:
            object = "s_%s*" % (self.name_sym_obj)
        if show_name == "surface":
            cmd.set("transparency", surf_transparency, object)
        cmd.disable("_sym_%s_%1.1f" % (self.name_sym_obj, self._cutoff))
    def get_model(self, selection=""):
        '''Return Model of selection.
        '''
        selection = _remove_brakets(selection)
        try:
            mol = cmd.get_model(selection)
        except:
            print("# Please specify a selection.")
            return 0
        return mol
class HP:
    def __init__(self):
        self.one_letter = {"ALA": ["A"],
                           "ARG": ["R"],
                           "ASN": ["N"],
                           "ASP": ["D"],
                           "ASX": ["B"],
                           "CYS": ["C"],
                           "CYH": ["C"],#protonated cys
                           "CYX": ["C"],#cystein
                           "GLN": ["Q"],
                           "GLU": ["E"],
                           "GLY": ["G"],
                           "GLX": ["Z"],
                           "HIS": ["H"],
                           "HIP": ["H"],#protonated
                           "HID": ["H"],#h on delta
                           "HIE": ["H"],#h on epsilon
                           "ILE": ["I"],
                           "LEU": ["L"],
                           "LYS": ["K"],
                           "MET": ["M"],
                           "MSE": ["M"],#Seleno- Methionin
                           "PHE": ["F"],
                           "PRO": ["P"],
                           "SER": ["S"],
                           "THR": ["T"],
                           "TRP": ["W"],
                           "TYR": ["Y"],
                           "VAL": ["V"],
                           "UNK": ["X"]}
        self.three_letter = {"A": ["ALA"],
                             "R": ["ARG"],
                             "N": ["ASN"],
                             "D": ["ASP"],
                             "B": ["ASX"],
                             "C": ["CYS", "CYH", "CYX"],
                             "Q": ["GLN"],
                             "E": ["GLU"],
                             "G": ["GLY"],
                             "Z": ["GLX"],
                             "H": ["HIS", "HIE", "HIP", "HID"],
                             "I": ["ILE"],
                             "L": ["LEU"],
                             "K": ["LYS"],
                             "M": ["MET"],
                             "MSE": ["MSE"],#selenomethionine
                             "F": ["PHE"],
                             "P": ["PRO"],
                             "S": ["SER"],
                             "T": ["THR"],
                             "W": ["TRP"],
                             "Y": ["TYR"],
                             "V": ["VAL"],
                             "X": ["UNK"]}
    def one2three(self, aa_list):
        '''Get one letter code from AA list.
        '''
        u = []
        tmp_list = [self.three_letter.get(l, ["UNK"]) for l in aa_list]
        for item in tmp_list:
            u.extend(item)
        return u
    def three2one(self, aa_list):
        '''Get one letter code from AA list.
        '''
        u = []
        tmp_list = [self.one_letter.get(l, ["X"]) for l in aa_list]
        for item in tmp_list:
            u.extend(item)
        return u
    def get_HP_set(self, hp_value="GES", output="names"):
        '''Get one letter code for hydrophobic,neutral,hydrophilic.
        USAGE:  get_HP_set(hp_value="GES")
                        hp_value: "GES"
                                  "KandD"
                                  "Rose"
                        output:   "names" :return AA names and category
                                  "values":return dic of AA names and values
        RETURNS: output names: (hydrophobic_AA,neutral_AA,hydrophilic_AA)
                 output values {'AA',value}
        '''
        hpvalues = {"GES": {"A": 1.6,
                            "C": 2,
                            "D": -9.2,
                            "E": -8.2,
                            "F": 3.7,
                            "G": 1,
                            "H": -3,
                            "I": 3.1,
                            "K": -8.8,
                            "L": 2.8,
                            "M": 3.4,
                            "N": -4.8,
                            "P": -0.2,
                            "Q": -4.1,
                            "R": -12.3,
                            "S": 0.6,
                            "T": 1.2,
                            "V": 2.6,
                            "W": 1.9,
                            "Y": -0.7, },
                    "KandD": {"A": 1.8,
                              "C": 2.5,
                              "D": -3.5,
                              "E": -3.5,
                              "F": 2.8,
                              "G": -0.4,
                              "H": -3.2,
                              "I": 4.5,
                              "K": -3.9,
                              "L": 3.8,
                              "M": 1.9,
                              "N": -3.5,
                              "P": -1.6,
                              "Q": -3.5,
                              "R": -4.5,
                              "S": -0.8,
                              "T": -0.7,
                              "V": 4.2,
                              "W": -0.9,
                              "Y": -1.3, },
                    "Rose": {"A": 0.74,
                             "C": 0.91,
                             "D": 0.62,
                             "E": 0.62,
                             "F": 0.88,
                             "G": 0.72,
                             "H": 0.78,
                             "I": 0.88,
                             "K": 0.52,
                             "L": 0.85,
                             "M": 0.85,
                             "N": 0.63,
                             "P": 0.64,
                             "Q": 0.62,
                             "R": 0.64,
                             "S": 0.66,
                             "T": 0.7,
                             "V": 0.86,
                             "W": 0.85,
                             "Y": 0.76, }}
        hydrophobic = {"GES":
                       ["F", "M", "I", "L", "V", "C", "W",
                           "A", "T", "G", "S", "MSE"],
                       "KandD":
                       ["A", "C", "I", "L", "M", "F", "V", "MSE"],
                       "Rose":
                       ["C", "I", "L", "M", "F", "W", "V", "MSE"]
                       }
        hydrophilic = {"GES": ["E", "K", "D", "R"],
                       "KandD":
                       ["R", "N", "D", "Q", "E", "H", "K"],
                       "Rose":
                       ["R", "N", "D", "Q", "E", "K", "P", "S"]
                       }
        neutral = {"GES": ["P", "Y", "H", "Q", "N"],
                   "KandD":
                       ["G", "P", "S", "T", "W", "Y"],
                       "Rose":
                       ["A", "G", "H", "T", "Y"]
                   }
        if output == "names":
            return (self.one2three(hydrophobic.get(hp_value)), self.one2three(neutral.get(hp_value)), self.one2three(hydrophilic.get(hp_value)))
        elif output == "values":
            return hpvalues.get(hp_values)
class set_views:
    def __init__(self):
        self.selection = "all"
    def mainchain(self, sele="all", obj_name="MC"):
        '''Create mainchain track in a Mainchain Object.
        mainchain(sele="all")
        '''
        group_name = obj_name + "_track_g"
        group_name = check_names(group_name)
        obj_name = obj_name + "_track"
        obj_name = check_names(obj_name)
        dist1 = check_names(obj_name + "_mcp_cont")
        dist2 = check_names(obj_name + "_scp_cont")
        self.selection = _remove_brakets(sele)
        try:
            cmd.create(obj_name, "(%s) and (n. C,N,O,CA)" % (self.selection))
            # create polar contacts of mainchain with preset and rename the
            # contacts
            cmd.dist(dist1, "(%s) and not (solvent or (polymer and not name n,o,h))" % self.selection,
                     "(%s) and not (solvent or (polymer and not name n,o,h))" % self.selection, quiet=1, mode=2, label=0, reset=1)
            cmd.enable(dist1)
            try:
                cmd.color(8, dist1)#color magenta
            except:
                pass
            cmd.show_as("sticks", obj_name)
            cmd.set("stick_radius", 0.35, obj_name)
            cmd.set("stick_transparency", 0.5, obj_name)
            cmd.dist(dist2, "(%s)" % self.selection, "(%s) and not (polymer and name n,o,h)" %
                     self.selection, quiet=1, mode=2, label=0, reset=1)
            cmd.enable(dist2)
            try:
                cmd.group(group_name, "%s %s %s" % (obj_name, dist1, dist2))
            except:
                pass
        except:
            print("# Can not create track. (no selection ?)")
        cmd.color("white", obj_name)
    def sidechain(self, sele="all"):
        '''Create sidechain.
        '''
        self.selection = _remove_brakets(sele)
        cmd.select("Sidechain", "(%s) and not n. C,N.O" % self.selection)
        cmd.show("lines", "Sidechain")
    def surface_atoms_preset(self, sele="all"):
        '''Get different important atoms.
        '''
        self.selection = _remove_brakets(sele, brakets=0)
        try:
            cmd.set("surface_color", -1, self.selection)
            cmd.select("HIS", "(%s) and (r. HIS,HIE,HID,HIP)" % self.selection)
            cmd.select(
                "HIS_ND1NE2", "(%s) and (r. HIS,HIE,HID,HIP) and (n. ND1,NE2)" % self.selection)
            cmd.color("gray", "(%s)" % self.selection)
            try:
                util.cba(29, "HIS_ND1NE2")
            except:
                pass
            cmd.select("ASP", "(%s) and (r. ASP)" % self.selection)
            cmd.select(
                "ASP_OD12", "(%s) and (r. ASP) and (n. OD1,OD2)" % self.selection)
            try:
                util.cba(29, "ASP_OD12")
            except:
                pass
            cmd.select("GLU", "(%s) and (r. GLU)" % self.selection)
            cmd.select(
                "GLU_OE12", "(%s) and (r. GLU) and (n. OE1,OE2)" % self.selection)
            try:
                util.cba(29, "GLU_OE12")
            except:
                pass
            cmd.select("ARG", "(%s) and (r. ARG)" % self.selection)
            cmd.select(
                "ARG_NH12", "(%s) and (r. ARG) and (n. NH1,NH2)" % self.selection)
            try:
                util.cba(29, "ARG_NH12")
            except:
                pass
            cmd.select("LYS", "(%s) and (r. LYS)" % self.selection)
            cmd.select("LYS_NZ", "(%s) and (r. LYS) and (n. NZ)" %
                       self.selection)
            try:
                util.cba(29, "LYS_NZ")
            except:
                pass
            cmd.select("CYS", "(%s) and (r. CYS,CYH,CYX)" % self.selection)
            cmd.select(
                "CYS_SG", "(%s) and (r. CYS,CYH,CYX) and (n. SG)" % self.selection)
            try:
                util.cba(29, "CYS_SG")
            except:
                pass
        except:
            print("# Could not create preset. (No Structure enabled?)")
def new_views(preset="track_mainchain", selection='sele'):
    '''Create Preset Views.
    '''
    #look if object is molecule and not a previous track run
    #names=[i for i in cmd.get_names("objects",enabled_only=1) if cmd.get_type(i)=="object:molecule"]
    selection = _remove_brakets(selection)
    if _check_selection(selection) == 0:
        print("# No (sele) selection, please select at least one atom.")
        return
    selectionname = selection
    selection = '(byres (%s))' % selection 
    if preset == "track_mainchain":
        #create for all enabled objects
        obj_name = check_names(selectionname)
        #try:
        #        names.remove(obj_name)
        #except ValueError:
        #        pass
        #selection=_remove_brakets(' or '.join(names))
        v = set_views()
        v.mainchain(selection, obj_name=obj_name)
    if preset == "surf_rel_res":
        #selection=_remove_brakets(' or '.join(names))
        #selection=_remove_brakets(selection)
        v = set_views()
        v.surface_atoms_preset(selection)
def set_surface(mode="standard"):
    '''Set surface to standard values.
    USAGE: 
            set_surface(mode="standard")
                    mode: "standard"
                          "symmetry"
    '''
    if mode == "standard":
        transparency = 0.0
        surfmode = -1
        selection = "all"
        cmd.set("surface_color", surfmode, selection)
        cmd.set("transparency", transparency, selection)
    elif mode == "symmetry":
        surf_transparency = 0.5
        selection = "all"
        set_surf_color(color="blue", sele="all")
        set_surf_color(color="white", sele="s_*")
        cmd.set("transparency", surf_transparency, selection)
    else:
        pass
def inv_enabled(selection='all'):
    '''Invert all enabled objects.
    USAGE: 
            inv_enabled selection="obj" # delete all enabled objects
            inv_enabled "selections" # delete all enabled selections
            inv_enabled "all" # delete all enabled 
    '''
    #invert
    selection = _remove_brakets(selection)
    print("# invert")
    vers = __pymol_version__
    if vers >= 1.21 and selection == "obj":
        selection = "objects"
    set_enabled = set(cmd.get_names(selection, enabled_only=1))
    set_all = set(cmd.get_names(selection, enabled_only=0))
    enabled = set_enabled
    disabled = set_all - set_enabled
    #print "disabled",disabled
    #cmd.disable(item)
    #cmd.enable(item)
    try:
        cmd.disable(' and '.join(enabled))
    except IndexError:
        pass
    try:
        cmd.enable(' and '.join(disabled))
    except IndexError:
        pass
def del_enabled(selection="obj"):
    '''Delete all enabled objects.
    USAGE: 
            del_enabled selection="obj" # delete all enabled objects
            del_enabled "selections" # delete all enabled selections
            del_enabled "all" # delete all enabled 
    '''
    selection = _remove_brakets(selection)
    vers = cmd.get_version()[1]
    if vers >= 1.21 and selection == "obj":
        selection = "objects"
    for item in cmd.get_names(selection, enabled_only=1):
        cmd.delete(item)
def rename_object(old_name="sele", new_name="new"):
    '''Rename any object without wizard.
    (The wizard cannot rename ramps and other stuff.)
    '''
    old_name = _remove_brakets(old_name)
    new_name = _remove_brakets(new_name)
    cmd.set_name(old_name, new_name)
def set_surf_color(color="-1", sele="all"):
    '''Set surface color idependently from structure.
    color: "-1" : default
    color: white
            "1"
            "2"
            "3"
            .
    '''
    color = _remove_brakets(color)
    sele = _remove_brakets(sele)
    try:
        color = int(color)
    except ValueError:
        pass
    cmd.set("surface_color", color, sele)
def get_pdb_lines(selection="all", state="0", enabled="0"):
    '''Get pdb lines from selection.
    USAGE: get_pdb_lines(selection="all",state="0",enabled="0")
    PYMOL: get_pdb (select),state
                    selection:   object or selection
                    state:       state
                    enabled:     1: get only enabled items
                                 0: get all items
    '''
    selection = _remove_brakets(selection)
    if selection != "all":
        names = [selection]
    else:
        names = [i for i in cmd.get_names(selection, enabled_only=int(enabled)) if cmd.get_type(
            i) == "object:molecule" or cmd.get_type(i) == "selection"]
    cnt = sum([cmd.count_atoms(i) for i in names])
    print("REMARK ________-Cut-Here-_________")
    print("REMARK Created by BNI Pymol Plugin v.%s" % (__version__))
    print("REMARK %s %s ATOMS" % (','.join(names), cnt))
    print("REMARK BEGIN PDB LINES FOR : %s" % (','.join(names)))
    for item in names:
        print(get_pdb(selection=item, state=int(state)))
        print("REMARK END   PDB LINES FOR : %s" % item)
    print("REMARK ________-Cut-END-_________")
def get_pdb(selection="all", state="0"):
    '''Get pdb lines from selection.
    USAGE: get_pdb(selection="all",state=0)
                    selection:   object or selection
                    state:       state
    '''
    selection = _remove_brakets(selection)
    pdbline = cmd.get_pdbstr(selection=selection, state=int(state))
    return pdbline
def del_sym(select_name="s_*"):
    '''Delete created symmetry mates.
    '''
    cmd.delete(select_name)
def select_hp(hp_value="KandD", selection="(sele)"):
    '''Select HP and color HP Values.
    '''
    u = HP()
    selection = _remove_brakets(selection)
    try:
        cmd.count_atoms(selection)
    except:
        print("# Please select at least one residue. (sele)")
        return
    (hydrophobic, neutral, hydrophilic) = u.get_HP_set(hp_value=hp_value)
    print("HP_Values %s" % hp_value)
    print("HP_pos (Hydrophobic) brown")
    print('\t' + ', '.join(hydrophobic))
    print("HP_ntr (Neutral) green")
    print('\t' + ', '.join(neutral))
    print("HP_neg (Hydrophilic) blue")
    print('\t' + ', '.join(hydrophilic))
    HP_pos = check_names("HP_pos")
    HP_ntr = check_names("HP_ntr")
    HP_neg = check_names("HP_neg")
    cmd.select(HP_pos, '%s and r. %s' % (selection, ','.join(hydrophobic)))
    cmd.color("brown", "(%s)" % HP_pos)
    cmd.select(HP_ntr, '%s and r. %s' % (selection, ','.join(neutral)))
    cmd.color("green", "(%s)" % HP_ntr)
    cmd.select(HP_neg, '%s and r. %s' % (selection, ','.join(hydrophilic)))
    cmd.color("blue", "(%s)" % HP_neg)
    try:
        HP_group = check_names(hp_value)
        members = "%s %s %s" % (HP_pos, HP_ntr, HP_neg)
        cmd.group(HP_group, members)
    except AttributeError: #no group command
        pass
def cgo_dic():
    cgo_dic = {
        'POINTS': 0.0,
        'LINES': 1.0,
        'LINE_LOOP': 2.0,
        'LINE_STRIP': 3.0,
        'TRIANGLES': 4.0,
        'TRIANGLE_STRIP': 5.0,
        'TRIANGLE_FAN': 6.0,
        #QUADS              : 7.0
        #QUAD_STRIP         : 8.0
        #POLYGON            : 9.0
        'STOP': 0.0,
        'NULL': 1.0,
        'BEGIN': 2.0,
        'END': 3.0,
        'VERTEX': 4.0,
        'NORMAL': 5.0,
        'COLOR': 6.0,
        'SPHERE': 7.0,
        'TRIANGLE': 8.0,
        'CYLINDER': 9.0,
        'LINEWIDTH': 10.0,
        'WIDTHSCALE': 11.0,
        'ENABLE': 12.0,
        'DISABLE': 13.0,
        'SAUSAGE': 14.0,
        'CUSTOM_CYLINDER': 15.0,
        'SHAPE_VERTEX': 16.0,
        'SHAPE_COLOR': 17.0,
        'SHAPE_NORMAL': 18.0,
        'FONT': 19.0,
        'FONT_SCALE': 20.0,
        'FONT_VERTEX': 21.0,
        'FONT_AXES': 22.0,
        'CHAR': 23.0,
        'ALPHA': 25.0, #transparency 1.00-0.01
        'LIGHTING': float(0x0B50)}
    return cgo_dic
def _remove_brakets(string, brakets=1):
    '''Remove brakets and quotation marks from string.
    brackets: 1: remove also brackets "(" and ")"
    '''
    if type(string) is not str:
        return string
    string = string.replace("'", "")
    string = string.replace('"', "")
    if brakets:
        string = string.replace("(", "")
        string = string.replace(")", "")
    return string
def get_center_atom(selection='all', name="patom", enabled="0", onlysele="0"):
    '''Get pseudo center of selection.
    USAGE: get_center_atom(selection,[name="center",enabled="0",onlysele=0])
    PYMOL: bni_get_center_atom selection,[name=center,enabled=0,onlysele=0])
    '''
    selection = _remove_brakets(selection)
    if selection != "all":
        names = [selection]
        #enabled="0"
    else:
        if onlysele == "0":
            names = [i for i in cmd.get_names(selection, enabled_only=int(enabled)) if cmd.get_type(
                i) == "object:molecule" or cmd.get_type(i) == "selection"]
        elif onlysele == "1":
            names = [i for i in cmd.get_names(
                selection, enabled_only=int(enabled)) if cmd.get_type(i) == "selection"]
        cnt = sum([cmd.count_atoms(i) for i in names])
    py_name = _remove_brakets(name)
    py_name = check_names(py_name)
    for item in names:
        print(item)
        center = getCenter(source=item)
        if center:
            if not name:
                py_name = "c_" + item
                py_name = check_names(py_name)
            get_atom(name=py_name, coord=center)
        else:
            print("# No Center to create.")
def getCenterCoord(coordlist):
    '''Get center coordinates of a coordinate list.
       coordlist=[(x,y,z),(x,y,z), ...]
    '''
    center = [0., 0., 0.]
    num_atom = 0
    for a in coordlist:
        num_atom = num_atom + 1
        center[0] = center[0] + a[0]
        center[1] = center[1] + a[1]
        center[2] = center[2] + a[2]
    if num_atom > 0:
        center[0] = center[0] / num_atom
        center[1] = center[1] / num_atom
        center[2] = center[2] / num_atom
    else:
        return
    return (center[0], center[1], center[2])
def getCenter(source='', logprint=1):
    '''
    DESCRIPTION:
        calculate the coordinate center of a selection
    USAGE:
        get_center selection
    '''
    source = _remove_brakets(source)
    try:
        mol = cmd.get_model(source)
    except:
        if logprint == 1:
            print("# Empty selection.")
            print("# create pseudo atom at (0.0,0.0,0.0)")
        return (0.0, 0.0, 0.0)
        #print getCenter.__doc__
    else:
        center = [0., 0., 0.]
        num_atom = 0
        for a in mol.atom:
            num_atom = num_atom + 1
            center[0] = center[0] + a.coord[0]
            center[1] = center[1] + a.coord[1]
            center[2] = center[2] + a.coord[2]
        if num_atom > 0:
            center[0] = center[0] / num_atom
            center[1] = center[1] / num_atom
            center[2] = center[2] / num_atom
        else:
            if logprint == 1:
                print("# No Atoms in selection.")
            return
        if logprint == 1:
            print("# gridcenter %8.3f %8.3f %8.3f\n" %
                  (center[0], center[1], center[2]))
        return (center[0], center[1], center[2])
def get_atom(name="atom", coord="(0.0,0.0,0.0)"):
    '''Create dummy atom on specified coordinate.
    USAGE: get_atom(name="atom",coord="(0.0,0.0,0.0)")
    PYMOL: get_atom name="atom",coord=(0.0,0.0,0.0)
    '''
    get_point(name, coord, radius="0.5", atom="1")
def check_names(name=''):
    '''Check if name is already in pymol viewer and create new name.
     USAGE: check_names(name='')
     OUTPUT: new name string: name_01
    '''
    names = cmd.get_names("all")
    if not names:
        return name
    check = 0
    for item in names:
        if item.lower() == name.lower():
            check = 1
            break
    if check == 1:
        try:
            a = cmd.get_unused_name(name + "_")
            return a
        except AttributeError:
            names = cmd.get_names("all")
            if not names:
                return name
            i = 1
            name_check = name
            check_again = 1
            a = 0
            while check_again:
                if check_again == 0:
                    break
                for item in names:
                    if item.lower() == name_check.lower():
                        try:#add additional number
                            a += 1
                            counter = int(name.split("_")[-1])
                            counter += a
                            name_check = "_".join(
                                name.split("_")[:-1]) + "_" + str(counter).zfill(2)
                            check_again = 1
                            i = 0 
                            break
                        except (IndexError, ValueError):#add new number
                            i += 1
                            name_check = name + "_" + str(i).zfill(2)
                            a = 0 # set count for additional number to 0
                            check_again = 1
                            break
                    else:
                        check_again = 0
            return name_check
    else:
        return name
def get_point(name="point", coord="(0.0,0.0,0.0)", radius="0.5", atom="0", state="new"):
    '''Create CGO Sphere point.
    USAGE: get_point(name="point",
                     coord="(0.0,0.0,0.0)",
                     radius="0.5",
                     atom="0")
    PYMOL: get_point name="point",             #Name of object
                     coord="(0.0,0.0,0.0)",    #coord of object
                     radius="0.5",             #radius of cgo
                     atom="0"                  #1: create dummy atom instead of cgo
                     state="new"               #"new": create unique name if name already exists
    '''
    radius = float(radius)
    if state == "new":
        name = check_names(name=name)
    if type(coord) is str:
        coord = _remove_brakets(coord)
        coord = coord.split(",")
        coord = tuple(map(float, coord))
        print(coord)
    else:
        print(coord)
    if int(atom) == 0:
        dic = cgo_dic()
        BEGIN = dic['BEGIN']
        TRIANGLES = dic['TRIANGLES']
        COLOR = dic['COLOR']
        SPHERE = dic['SPHERE']
        LINES = dic['LINES']
        NORMAL = dic['NORMAL']
        VERTEX = dic['VERTEX']
        END = dic['END']
        p = [SPHERE,
             coord[0], coord[1], coord[2], radius]
        cmd.load_cgo(p, name)
    elif int(atom) == 1:
        record_name = 'ATOM  '
        atom_id = 5
        atom_name = ' PO '
        alternate = ' '
        residue_name = "INT"
        chain_id = 'P'
        residue_id = 1
        icode = " "
        x = coord[0]
        y = coord[1]
        z = coord[2]
        occ = 0.00
        bfac = 100.00
        free = 6 * ' '
        seqid = "POIN"
        symbol = " T"
        charge = " "
        state = -1
        cmd.read_pdbstr('%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%6s%4s %2s%2s' %
                        (record_name,
                         atom_id,
                         atom_name,
                         alternate,
                         residue_name,
                         chain_id,
                         residue_id,
                         icode,
                         x,
                         y,
                         z,
                         occ,
                         bfac,
                         free,
                         seqid,
                         symbol,
                         charge
                         ), name, state)
        try:
            cmd.show_as("nb_spheres", name)
        except AttributeError:
            print(
                "# cmd.as() has changed to cmd.show_as() in later and previous versions of PyMOL ")
            #cmd.show_as("nb_spheres",name)
def get_axes(name="axes", origin="(0.0,0.0,0.0)", lenght="10.0", linewidth="3.0"):
    '''Create Axes.
    USAGE: get_axes(name="axes",origin="(0.0,0.0,0.0)",lenght="10.0",linewidth="1.0")
    PYMOL: axes name="axes",origin="(0.0,0.0,0.0)",lenght="10.0",linewidth="1.0"
                name:  name of object
                origin:  zero coordinates of axes or "box" for origin coordinates of selection min coordinates
                lenght:  lenght of axes
                linewidth: thickness of axes
    '''
    if origin == "box":
        origin, (xmax, ymax, zmax) = cmd.get_extent(selection)
    try:
        linewidth = float(linewidth)
    except ValueError:
        print("# Wrong linewidth input 1.0 is used instead.")
        linewidth = 1.0
    if type(origin) is str:
        origin = _remove_brakets(origin)
        origin = origin.split(",")
        origin = tuple(map(float, origin))
    else:
        pass
    lenght = float(lenght)
    dic = cgo_dic()
    BEGIN = dic['BEGIN']
    TRIANGLES = dic['TRIANGLES']
    COLOR = dic['COLOR']
    SPHERE = dic['SPHERE']
    LINES = dic['LINES']
    NORMAL = dic['NORMAL']
    VERTEX = dic['VERTEX']
    END = dic['END']
    LINEWIDTH = dic['LINEWIDTH']
    a = [LINEWIDTH, linewidth,
         BEGIN,
         LINES,
         COLOR, 0.0, 1.0, 0.0,
         VERTEX,
         0.0 + origin[0], 0.0 + origin[1], 0.0 + origin[2],
         VERTEX,
         0.0 + origin[0], lenght + origin[1], 0.0 + origin[2],
         COLOR, 1.0, 0.0, 0.0,
         VERTEX,
         0.0 + origin[0], 0.0 + origin[1], 0.0 + origin[2],
         VERTEX,
         lenght + origin[0], 0.0 + origin[1], 0.0 + origin[2],
         COLOR, 0.0, 0.0, 1.0,
         VERTEX,
         0.0 + origin[0], 0.0 + origin[1], 0.0 + origin[2],
         VERTEX,
         0.0 + origin[0], 0.0 + origin[1], lenght + origin[2],
         END,
         COLOR, 1.0, 0.0, 0.0,
         SPHERE,
         0.0 + origin[0], 0.0 + origin[1], 0.0 + origin[2], 0.2]
    print("# Axes\n-----\n\tX :red\n\tY :green\n\tZ :blue")
    cmd.load_cgo(a, name)
def get_box_handler(selection="(sele)", app=None, action="medium", name="box", cushion="1.5", draw="full", coord="None", linewidth="3.0"):
    '''Get a box using get_box GUI.
    '''
    action = _remove_brakets(action)
    name = check_names(name)
    cushion = cushion
    draw = str(draw)
    linewidth = str(linewidth)
    if action == "set" and app:
        value = tkSimpleDialog.askstring('Box settings',
                                         'Specify [cushion] [draw full/lines] [linewidth]',
                                         parent=app.root)
        value = _remove_brakets(value)
        if value:
            value = value.split()
        else:
            return
        if len(value) == 1:
            cushion = value[0]
        elif len(value) == 2:# width / dpi
            cushion = value[0]
            draw = str(value[1])
        elif len(value) >= 3:#width width_unit
            cushion = value[0]
            draw = str(value[1])
            linewidth = str(value[2])
        else:
            print("# Wrong input. Set at least size.")
            return
        try:
            float(cushion)
        except ValueError:
            print("# Cushion has to be a number.")
            return
        units = ["full", "lines"]
        if not draw.lower().strip() in units:
            print("# Draw has to be either %s " % (" or ".join(units)))
            return
        try:
            float(linewidth)
        except ValueError:
            print("# Linewidth has to be a number.")
            return
    elif action == "small":
        cushion = 1
        draw = "full"
        linewidth = 3
    elif action == "medium":
        cushion = 2
        draw = "full"
        linewidth = 3
    elif action == "large":
        cushion = 5
        draw = "full"
        linewidth = 3
    get_box(selection=selection, name=name,
            cushion=cushion, draw=draw, linewidth=linewidth)
def bni_group(group_name="TEST_Group"):
    '''Check if group is present if so add new group.
    '''
    #try:
    name = check_names(group_name)
    try:
        cmd.group(name)
    except:
        return
    return name 
def get_box(name="box", selection="all", cushion="1.5", draw="full", coord="None", linewidth="3.0"):
    '''Create a box.
    USAGE: get_box(name="box",
                      selection="all",
                      cushion=0.0,
                      draw="full",
                      coord="None" or "(xmin,ymin,zmin,xmax,ymax,zmax)",
                      linewidth="1.0"
                      )
    PYMOL: get_box name=box,selection=all,cushion=2.0,draw=full,coord=None,linewidth=1.0
              name:       name of the cgo box
              selection:  selection for min and max of box calculation
              cushion  :  extend box with cushion in Angstroms each side
                          or each side separately by "(xcmin,ycmin,zcmin,xcmax,ycmax,zcmax)"
              draw:     lines:       draw a box with lines
                        full:        draw a full box with each side separated
                                     to hide box walls separately
              coord    :  if coord is given as "(xmin,ymin,zmin,xmax,ymax,zmax)"
                          selection will be ignored
              linewidth: 1.0 : in case of box with lines define linewidth
    '''
    draw = _remove_brakets(draw)
    name = _remove_brakets(name)
    gn = bni_group("BNI-box")
    dic = cgo_dic()
    try:
        cushion = float(cushion)
        cushion = (cushion, cushion, cushion, cushion, cushion, cushion)
    except ValueError:
        try:
            cushion = _remove_brakets(cushion)
            cushion = cushion.split(",")
            cushion = tuple(map(float, cushion))
        except ValueError:
            pass
        if len(cushion) == 6:
            pass
        else:
            cushion = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    linewidth = float(linewidth)
    if coord == "None":
        (xmin, ymin, zmin), (xmax, ymax, zmax) = cmd.get_extent(selection)
        (xmin, ymin, zmin), (xmax, ymax, zmax) = (xmin - cushion[0], ymin - cushion[
            1], zmin - cushion[2]), (xmax + cushion[3], ymax + cushion[4], zmax + cushion[5])
    else:
        if type(coord) is str:
            coord = _remove_brakets(coord)
            coord = coord.split(",")
            coord = tuple(map(float, coord))
        else:
            pass
        xmin, ymin, zmin, xmax, ymax, zmax = coord
        xmin, ymin, zmin, xmax, ymax, zmax = xmin - cushion[0], ymin - cushion[
            1], zmin - cushion[2], xmax + cushion[3], ymax + cushion[4], zmax + cushion[5]
    edge01 = (xmax, ymax, zmax)
    edge02 = (xmin, ymax, zmax)
    edge03 = (xmin, ymax, zmin)
    edge04 = (xmax, ymax, zmin)
    edge05 = (xmax, ymin, zmin)
    edge06 = (xmax, ymin, zmax)
    edge07 = (xmin, ymin, zmax)
    edge08 = (xmin, ymin, zmin)
    normx = (-1.0, 0.0, 0.0)
    normx1 = (1.0, 0.0, 0.0)
    normy = (0.0, -1.0, 0.0)
    normy1 = (0.0, 1.0, 0.0)
    normz = (0.0, 0.0, -1.0)
    normz1 = (0.0, 0.0, 1.0)
    if draw == "lines":
        BEGIN = dic['BEGIN']
        LINES = dic['LINES']
        LINEWIDTH = dic['LINEWIDTH']
        VERTEX = dic['VERTEX']
        END = dic['END']
        box = [
            LINEWIDTH, linewidth,
            BEGIN,
            LINES,
            VERTEX, edge01[0], edge01[1], edge01[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            VERTEX, edge01[0], edge01[1], edge01[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            VERTEX, edge01[0], edge01[1], edge01[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            END
        ]
        name = check_names(name)
        cmd.load_cgo(box, name)
        try:
                        #name=name+"_%d"%(i+1)
                        #name=check_names(name)
            cmd.group(gn, name)
        except:
            pass
    elif draw == "full":
        box = []
        BEGIN = dic['BEGIN']
        TRIANGLES = dic['TRIANGLES']
        COLOR = dic['COLOR']
        SPHERE = dic['SPHERE']
        LINES = dic['LINES']
        NORMAL = dic['NORMAL']
        VERTEX = dic['VERTEX']
        END = dic['END']
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge01[0], edge01[1], edge01[2],
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            NORMAL, normy[0], normy[1], normy[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            END
        ]
        box.append(a)
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge01[0], edge01[1], edge01[2],
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            NORMAL, normx[0], normx[1], normx[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            END
        ]
        box.append(a)
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            NORMAL, normz1[0], normz1[1], normz1[2],
            VERTEX, edge04[0], edge04[1], edge04[2],
            END
        ]
        box.append(a)
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge03[0], edge03[1], edge03[2],
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            NORMAL, normx1[0], normx1[1], normx1[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            END
        ]
        box.append(a)
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge01[0], edge01[1], edge01[2],
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge02[0], edge02[1], edge02[2],
            NORMAL, normz[0], normz[1], normz[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            END
        ]
        box.append(a)
        a = [
            BEGIN, TRIANGLES,
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge06[0], edge06[1], edge06[2],
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge05[0], edge05[1], edge05[2],
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge08[0], edge08[1], edge08[2],
            NORMAL, normy1[0], normy1[1], normy1[2],
            VERTEX, edge07[0], edge07[1], edge07[2],
            END
        ]
        box.append(a)
        color = ["wheat",
                 "lightpink",
                 "palecyan",
                 "darksalmon",
                 "gray90",
                 "palegreen"]
        for i, (item, color) in enumerate(zip(box, color)):
            n = name + "-%s" % str(i + 1).zfill(2)
            n = check_names(n)
            cmd.load_cgo(item, n)
            cmd.color(color, n)
            try:
                cmd.group(gn, n)
            except:
                pass
                #print("Cannot group %s"%group_name)
def fast_sym(selection='', cutoff=4.0, delete_sym=1, use=1, app=None):
    '''Create a view of symmetry involved interactions.
    USAGE: fast_sym selection,cutoff=4.0,delete_sym=1
                    selection:   object for symmetry calculation
                    cutoff:      distance cutoff
                    delete_sym:  delete previous sym files ("s_*")
    '''
    selection = _remove_brakets(selection)
    if cutoff == "set" and app:
        cutoff = tkSimpleDialog.askstring('Input Symetry Distance',
                                          'Please specify a Distance e.g 4.25',
                                          parent=app.root)
    if cutoff:
        try:
            cutoff = float(cutoff)
            use = int(use)
            delete_sym = int(delete_sym)
        except:
            print("# Wrong input")
            print(fast_sym.__doc__)
        else:
            cryst = Sym()
            if delete_sym == 1:#delete previous symetry objects
                print(delete_sym)
                del_sym()
            if selection:
                if cmd.get_object_list(selection) is None:
                    print(
                        "# No structure in selection (sele). Select at least one atom")
                    print("# in a structure.")
                    return
                if len(cmd.get_object_list(selection)) > 1:
                    print("# WARNING: more than one object in selection.")
                    print(
                        "# WARNING: first object will be used for calculation")
                try:
                    selection = selection + " and " + \
                        cmd.get_object_list(selection)[0]#just use first hit
                except IndexError:
                    print("# no object enabled")
                    print("# enable an object")
                    selection = ''
                    return
            elif selection == '':#get first object
                names = cmd.get_names(enabled_only=1)
                #selection=' or '.join(names)
                if len(names) > 1:
                    print("# WARNING: more than one object is enabled")
                    print("# WARNING: first object is used")
                elif len(names) == 0:
                    print("# load a structure into PyMOL first.")
                    return
                try:
                    selection = names[0]#just use first hit
                except IndexError:
                    print("# no object enabled")
                    print("# enable an object")
                    selection = ''
            else:
                print("# load a structure into pymol. No structure selected")
                return
            mol = cryst.get_model(selection)
            if mol:
                if len(mol.atom) == 0:
                    print("# no atoms in selection")
                elif use == 1:
                    original = selection
                    cmd.hide("surface")
                    sym_select = cryst.create_sym(
                        cutoff=cutoff, selection=selection)
                    set_surface(mode="symmetry")
                    org_select = cryst.create_org(cutoff=cutoff, sele="s_*")
                    cryst.show_org(
                        select_name=org_select, show_name="surface", cutoff=cutoff, surf_transparency=0.5)
                    cryst.show_org(
                        select_name=org_select, show_name="stick", cutoff=cutoff, surf_transparency=0.5)
                    cryst.show_sym(select_name=sym_select, show_name="surface",
                                   cutoff=cutoff, surf_transparency=0.5, object="s_*")
                    cryst.show_sym(select_name=sym_select, show_name="stick",
                                   cutoff=cutoff, surf_transparency=0.5, object="s_*")
                cmd.delete("tmp_org")

import Pmw
class multi_pdb:
    '''Multiple input load. 
    '''
    filedir = '.'
    def __init__(self, app, case=1, message=None):
        self._message = message
        try:
            self._app = app.root
        except:
            self._app = app
        self.filedir = self.__class__.filedir
        self.filename = []
        self.filetypes = [('PDB File', '.pdb'),
                          ('PDB File', '.ent'),
                          ('Modeller PDB File', '.B9999*.pdb'),
                          ('Modeller PDB File', '.B9999*_fit.pdb'),
                          ('Delphi PHI File', '.phi'),
                          ('all', '*'),
                          ]
    def _map(self, filenamelist):
        '''load Electon density map.
        '''
        pass
    def _casox_map(self, filenamelist, repr="mesh"):
        '''Load a CaSoX map (xplor format) and generate different level views.
         USAGE: _casox_map(filenamelist,repr="mesh")
                     repr:   "mesh"   :show levels in mesh format
                             "surf"   :show levels in surf format
         ADDITIONAL: this will create a map and different mesh representations
                     in different colors depending on casox ligsite values
                     level 0 --> "white"       --> casox ligsite level 1
                     level 1 --> "bluewhite"   --> casox ligsite level 2
                     level 2 --> "palecyan"    --> casox ligsite level 3
                     level 3 --> "lightteal"   --> casox ligsite level 4
                     level 4 --> "lightorange" --> casox ligsite level 5
                     level 5 --> "orange"      --> casox ligsite level 6
                     level 6 --> "red"         --> casox ligsite level 7 (complete closed cavity)
                     level-1 -->  "blue"       --> casox protein level 0 grid hull                   
        '''
        if len(filenamelist) != 1:
            print("# please select just one CASoX map.")
            file_path = norm_path(filenamelist[0])
            self.__class__.filedir = file_path.dirname
            return
        else:
            file_path = norm_path(filenamelist[0])
            objcasox = file_path.splitname + "_map"
            self.open_filename(file_path.all, obj=objcasox, format="xplor")
            self.__class__.filedir = file_path.dirname
        if cmd.get_type(objcasox) == "object:map":
            levels = [-1, 6, 5, 4, 3, 2, 1, 0]
            colors = ["blue", "red", "orange", "lightorange",
                      "lightteal", "palecyan", "bluewhite", "white"]
            print("# CaSoX Map")
            print("# --------------")
            print("# Loaded by BNI-Tools v.%s" % __version__)
            print("# Map: %s " % (objcasox))
            print("# Levels for Ligsite 1 - 7(closed cavity)")
            cmd.feedback("disable", "cmd", "everything")
            try:
                for level, color in zip(levels, colors):
                    objmesh = file_path.splitname + "_%d" % (level + 1)
                    if repr == "mesh":
                        cmd.isomesh(objmesh, objcasox, 1.0)
                    elif repr == "surf":
                        cmd.isosurface(objmesh, objcasox, 1.0)
                        cmd.hide("surface", objmesh)
                    cmd.isolevel(objmesh, float(level))
                    cmd.color(color, objmesh)
                    print("%15s : %15s" % (objmesh, color))
            except:
                print("# casox map could not be recognized")
                print("# loading of casox map failed.")
            print("# --------------")
            cmd.feedback("enable", "cmd", "everything")
    def _delphi(self, filenamelist):
        '''Load phi and pdb file to create a delphi surface.
        '''
        if len(filenamelist) > 2:
            print("# Select .phi or .dx map and corresponding .pdb file.")
            file_path = norm_path(filenamelist[0])
            self.__class__.filedir = file_path.dirname
            return
        else:
            objphi = "e_lvl"
            map = 0
            for item in filenamelist:
                file_path = norm_path(item)
                if file_path.splitext.lower() == ".phi" or file_path.splitext.lower() == ".dx":
                    objmap = file_path.splitname + "_map"
                    objmapext = file_path.splitext.lower()
                    objphi = file_path.splitname
                    self.open_filename(item, obj=objmap)
                    map = 1
                else:
                    objpdb = file_path.splitname
                    self.open_filename(item, obj=objpdb)
            self.__class__.filedir = file_path.dirname
            if map == 1 and cmd.get_type(objmap) == "object:map":
                objramp = objphi + "_ramp"
                cmd.ramp_new(objramp, objmap, (-15, 0, +15))
                if objmapext == ".phi":
                    print("# DELPHI PHI Map")
                elif objmapext == ".dx":
                    print("# DX Map")
                print("# --------------")
                print("# Loaded by BNI-Tools v.%s" % __version__)
                print("# Map    : %s" % objmap)
                print("# Ramp   : %s" % objramp)
                if len(filenamelist) == 2:
                    cmd.set("surface_color", objramp, objpdb)
                    print("# File   : %s" % objpdb)
                else:
                    print("# load pdb file as object into pymol and")
                    print("# set surface_color,%s,(object)" % objramp)
                print("#--------------")
                print(
                    "# to change ramp values use str+mid mouse button on ramp")
                print("# or use")
                print("# ramp_new %s,%s,(min,mid,max)" % (objramp, objmap))
                print("# to set ramp values.")
            else:
                print("# No Map loaded ... (extention has to be .phi or .dx)")
    def _modeller(self, filenamelist):
        '''Load multiple modeller files into Pymol.
        '''
        filenamelist = list(filenamelist)
        filenamelist.sort()
        filenamelist = tuple(filenamelist)
        self.filename = []
        for item in filenamelist:
            file_path = norm_path(item)
            self.filename.append(file_path.basename)
        self.__class__.filedir = file_path.dirname
        objfunc = []
        seqid = []
        modvers = []
        filelist = []
        object = []
        for files in self.filename:
            file_path = norm_path(self.__class__.filedir, files)
            filelist.append(file_path.all)
            obj = file_path.splitname
            if ".B9999" in obj: 
                obj = obj.split(".B9999")
                obj1 = obj[0]
                try:
                    obj2 = obj[1]
                except IndexError:
                    obj2 = ""
                    obj = obj1
                try:
                    obj2 = int(obj2)
                except ValueError:
                    pass
                obj = obj1 + '_' + "M" + str(obj2)
            elif ".BL" in obj:
                obj = obj.split(".BL")
                obj1 = obj[0]
                try:
                    obj2 = obj[1][:-4]
                except IndexError:
                    obj2 = ""
                    obj = obj1
                try:
                    obj2 = int(obj2)
                except ValueError:
                    pass
                obj = obj1 + '_' + 'L' + str(obj2)
            object.append(obj)
            try:
                fileobj = open(file_path.all, "r")
            except IOError:
                print("# cannot read modeller file  %s " % file_path.basename)
                object.pop() 
                continue
            else:
                sid = 0.0
                objf = 0.0
                while 1:
                    a = fileobj.readline()
                    if not a:
                        break
                    if "ATOM" in a[:4] or "HETATM" in a[:6]:
                        break
                    if "EXPDTA" in a[:6] and "MODELLER" in a:
                        continue
                    elif "REMARK" in a[:6] and "MODELLER OBJECTIVE FUNCTION" in a:
                        try:
                            objf = float(a.split()[-1])
                        except (IndexError, ValueError):
                            print(
                                "# Cannot parse objective funktion for %s" % file_path.basename)
                            pass
                    elif "REMARK" in a[:6] and "MODELLER BEST TEMPLATE" in a:
                        try:
                            sid = float(a.split()[-1])
                        except (IndexError, ValueError):
                            print("# cannot parse SEQ ID for %s" %
                                  file_path.basename)
                            pass
                objfunc.append(objf)
                seqid.append(sid)
        s = list(zip(objfunc, filelist, object, seqid))
        s.sort()
        print("# Modeller Files")
        print("# --------------")
        print("# Loaded by BNI-Tools v.%s" % __version__)
        print("# Sorted by Modeller Objective Function")
        print("# Name ;Objective Function; Seq-ID")
        for objfunc, item, obj, seqid in s:
            if "%.2f" % objfunc == "0.00":
                self.open_filename(item, obj=obj)
            else:
                self.open_filename(item, obj=obj, title="%.2f" % objfunc)
            print("%-15s ; %.2f ;  %.2f%%" % (obj, objfunc, seqid))
        print("# --------------")
    def _amber(self, filenamelist):
        '''Load multiple amber pdb files with info file and energy.
        '''
        filenamelist = list(filenamelist)
        filenamelist.sort()
        filenamelist = tuple(filenamelist)
        for item in filenamelist:
            file_path = norm_path(item)
            self.filename.append(file_path.basename)
        self.__class__.filedir = file_path.dirname
        filelist = []
        print("# Amber Files")
        print("# --------------")
        print("# Loaded by BNI-Tools v.%s" % __version__)
        info_file_ext = ".info"
        print("_" * 60)
        print("%-40s:%-20s" % ("Name", "Energy"))
        print("-" * 60)
        for files in self.filename:
            file_path = norm_path(self.__class__.filedir, files)
            filelist.append(file_path.all)
            obj = file_path.splitname
            try:
                self.open_filename(file_path.all, obj=obj)
            except IOError:
                print("# cannot open file %s" % file_path.basename)
                continue
            try:
                info_file = norm_path(
                    file_path.dirname, file_path.splitname + info_file_ext)
                if os.path.isfile(info_file.all):
                    fileobj = open(info_file.all, "r")
                #try to find info file if not same name as pdb file
                elif "_min." in file_path.basename:
                    info_name = file_path.splitname.split("_min")[0]
                    info_file = norm_path(
                        file_path.dirname, info_name + info_file_ext)
                    fileobj = open(info_file.all, "r")
                else:
                    print("# %s has no amber .info file in same directory." %
                          (file_path.basename))
                    print("# therefore no energy could be read in.")
                    continue
            except IOError:
                print("# cannot read amber file  %s " % file_path.basename)
                continue
            else:
                energy = 0.0
                energy_sign = "NSTEP"
                remark_text = []
                begin_read = 0
                begin_energy = 0
                while 1:
                    a = fileobj.readline()
                    if not a:
                        break
                    if begin_energy == 1:
                        try:
                            energy = float(a.split()[1])
                        except (IndexError, ValueError):
                            energy = 0.0
                        begin_energy = 0
                    if energy_sign in a:
                        begin_read = 1
                        begin_energy = 1
                        remark_text.append(a)
                fileobj.close()
                print("%-40s:%-20.3f" % (obj, energy))
                try:
                    cmd.set_title(obj, -1, "%-5.1f" % (energy))
                except:
                    pass
        print("-" * 60)
    def _autodock(self, filenamelist):
        '''Load multiple autodock MODEL files.
        '''
        filenamelist = list(filenamelist)
        filenamelist.sort()
        filenamelist = tuple(filenamelist)
        for item in filenamelist:
            file_path = norm_path(item)
            self.filename.append(file_path.basename)
        self.__class__.filedir = file_path.dirname
        filelist = []
        print("# AutoDock Files")
        print("# --------------")
        print("# Loaded by BNI-Tools v.%s" % __version__)
        print("_" * 60)
        print("%-30s:%-30s" % ("Name", "Nr. of Models"))
        print("-" * 60)
        for files in self.filename:
            file_path = norm_path(self.__class__.filedir, files)
            filelist.append(file_path.all)
            obj = file_path.splitname
            try:
                self.open_filename(file_path.all, obj=obj)
            except IOError:
                print("# cannot open file %s" % file_path.basename)
                continue
            try:
                fileobj = open(file_path.all, "r")
            except IOError:
                print("# cannot read autodock file  %s " % file_path.basename)
                continue
            else:
                cluster = []
                cluster_sign = "Number of conformations in this cluster"
                cluster_split = "="
                energy = []
                energy_sign = "Estimated Free Energy of Binding"
                energy_split = "="
                remarks = []
                remark_text = []
                model_state = []
                model_state_id = 0
                model_id = []
                begin_read = 0
                while 1:
                    a = fileobj.readline()
                    if not a:
                        break
                    if "ENDMDL" in a[:6] and begin_read == 1:
                        begin_read = 0
                        remarks.append(remark_text)
                        remark_text = []
                        continue
                    if "REMARK" in a[:6] or "USER" in a[:4] and begin_read == 1:
                        try:
                            remark_text.append(a.split()[1])
                        except IndexError:
                            remark_text.append(" ")
                        if cluster_sign in a:
                            try:
                                cluster.append(
                                    int(a.split(cluster_split)[1].split()[0]))
                            except IndexError:
                                cluster.append(0)
                        if energy_sign in a:
                            try:
                                energy.append(
                                    float(a.split(energy_split)[1].split()[0]))
                            except IndexError:
                                energy.append(0.0)
                        continue
                    if "MODEL" in a[:5]:
                        try:
                            model_id.append(int(a.split()[-1]))
                        except IndexError:
                            model_id.append(0)
                        model_state_id += 1
                        model_state.append(model_state_id)
                        begin_read = 1
                fileobj.close()
                print("%-26s:%-30d" % (obj, len(model_state)))
                print("-" * 60)
                print("%-5s%-4s:%-4s:%-10s :%-10s" %
                      (" ", "Nr", "Id", "Energy", "Cluster"))
                for item in zip(model_state, model_id, energy, cluster):
                    #
                    print("     %-4d:%-4d:%10.3f : %-3s %s" %
                          (item + ("*" * int(item[-1]),)))
                    try:
                        cmd.set_title(
                            obj, item[0], "%5.2f [%2d]" % (item[2], item[3]))
                    except:
                        pass
                print("-" * 60)
                if file_path.splitext == ".dlg":
                    print(
                        "# DLG File loaded (best docking mode in cluster is shown.)")
                    # additional read all docked models
        #s=zip(objfunc,filelist,object,seqid)
        #s.sort()
    def _oneobject(self, filenamelist, sorting=1):
        '''Load multiple files into one object.
           Object is named by first filename.
           If yasara table (.tbl) is present in selection try to read energies.
           sorting: 1: sort by filename
        '''
        filenamelist = list(filenamelist)
        if sorting == 1:
            filenamelist.sort()
        tablename = None
        try:
            table_index = [
                filenamelist.index(a) for a in filenamelist if "tab" in norm_path(a).splitext]
            if table_index:
                for tab_index in table_index:
                    tablename = filenamelist.pop(tab_index)
                tablename = norm_path(tablename)
        except IndexError:
            pass
        try:
            object_name = norm_path(filenamelist[0]).splitname
        except IndexError:
            object_name = "obj1"
        filenamelist = tuple(filenamelist)
        for item in filenamelist:
            file_path = norm_path(item)
            self.filename.append(file_path.basename)
        self.__class__.filedir = file_path.dirname
        for files in self.filename:
            file = norm_path(self.__class__.filedir, files).all
            self.open_filename(file, obj=object_name)
        if tablename:
            try:
                fileobj = open(tablename.all, "r")
            except IOError:
                print("# cannot read table file  %s " % tablename.basename)
            else:
                begin_read = 0
                energy = []
                while 1:
                    a = fileobj.readline()
                    if not a:
                        break
                    if begin_read == 1:
                        energy.append(a.strip().split())
                    if "snapshot" in a.strip().lower():
                        print("Sanapshot   Energy")
                        print("------------------")
                        begin_read = 1
                        continue
                fileobj.close()
                energy = [e for e in energy if len(e) == 2]
                for i, (snap, value) in enumerate(energy):
                    print("%-5s:%-30s" % (snap, value))
                    try:
                        cmd.set_title(
                            object_name, i + 1, "%s:%s" % (snap, value))
                    except:
                        pass
                print("-" * 60)
    def _all(self, filenamelist):
        '''Load multiple files.
        '''
        filenamelist = list(filenamelist)
        filenamelist.sort()
        filenamelist = tuple(filenamelist)
        for item in filenamelist:
            file_path = norm_path(item)
            self.filename.append(file_path.basename)
        self.__class__.filedir = file_path.dirname
        for files in self.filename:
            file = norm_path(self.__class__.filedir, files).all
            self.open_filename(file, obj="None")
    def fetchFileName(self, ext="None"):
        '''Get multiple filenmes loaded into Pymol with special filename handling.
        USAGE: fetchFileName(ext="None")
                 ext:   "None"     :load files with filename extention
                        "modeller" :load modeller files sorted
                                    by modeller objective function
                                    and rename filenames .B9999000X to _X
                        "autodock" :load autodock pdb files (docked MODELs)
                                    with energy and cluster in title
                        "ambermin" :load amber minimization pdb with energy
                                    in .info files
                        "delphi"   :load delphi map (*.phi) and pdb file for surface representation
                                    or abps map (*.dx) and pdb file for surface representation
                        "casox"    :load CaSoX Map (xplor like Soft Ligsite Cavity Explorer)
        '''
        title = 'Open File(s)'
        if ext == "None":
            filetypes = self.filetypes
        elif ext == "oneobject":
            filetypes = self.filetypes
        elif ext == "modeller":
            title = 'Select Modeller PDB file (.B9999*.pdb)'
            filetypes = [('Modeller PDB File', '.B9999*.pdb'), ('all', '*')]
        elif ext == "autodock":
            title = 'Select Autodock .dlg file or .dlg.pdb (old) file and receptor pdb.'
            filetypes = [('AutoDock DLG file', '.dlg'), ('AutoDock docked PDB', '.pdb'), (
                'AutoDock DLG file', '.pdb'), ('AutoDock docked PDB', '.dlg.pdb'), ('all', '*')]
        elif ext == "ambermin":
            title = 'Select amber minimization pdb file (.info file in same directory displays energies).'
            filetypes = [('Amber minimization', '.pdb'), ('all', '*')]
        elif ext == "delphi":
            title = 'Select map and also the corresponding .pdb file'
            filetypes = [('DelPhi PHI Map', '.phi'), ('DelPhi PHI Map', '.pdb'),
                         ('ABPS dx Map', '.pdb'), ('ABPS dx Map', '.dx'), ('all', '*')]
        elif ext == "casox":
            title = 'Select CASoX Map'
            filetypes = [
                ('CASoX Map', '.omap'), ('CASoX Map', '.map'), ('all', '*')]
        filenamelist = tkFileDialog.askopenfilenames(title=title,
                                                     initialdir=self.__class__.filedir,
                                                     filetypes=filetypes,
                                                     )
        self.filename = []
        self.filedir = self.__class__.filedir
        if __pymol_version__ > 1.3:
            try:
                filenamelist = self._app.splitlist(filenamelist)
            except:
                print(
                    "# cannot create workaround for python askopenfielanmes bug.")
        if filenamelist:
            if ext == "None":
                self._all(filenamelist)
            elif ext == "modeller":
                self._modeller(filenamelist)
            elif ext == "autodock":
                self._autodock(filenamelist)
            elif ext == "ambermin":
                self._amber(filenamelist)
            elif ext == "delphi":
                self._delphi(filenamelist)
            elif ext == "casox":
                self._casox_map(filenamelist)
            elif ext == "oneobject":
                self._oneobject(filenamelist)
    def open_filename(self, file, obj="None", state="None", format="None", finish="None", discrete="None", title="None", allinone="0"):
        '''Open file in pymol.
        USAGE: open_filename(file,
                             obj="None",
                             state="None",     # not implemented
                             format="None",    # specify format
                             finish="None",    # not implemented
                             discrete="None",  # not implemented
                             title="None"      # write title
                             allinone="0")     # all in one object (first object name) with different states 
        '''
        if obj == "None":
            obj = norm_path(file).splitname
        else:
            pass
        if state == "None":
            state = -1
        if format == "None":
            cmd.load(file, object=obj, state=state)
        else:
            cmd.load(file, object=obj, state=state, format=format)
            #print "loaded",file,obj,format,state
        if title == "None":
            title = ""
        else:
            cmd.set_title(obj, state, title)
def norm_path(pathname, *join, **abs):
    '''Create normalized path and join path items.
    USAGE: norm_path(pathname
                     [,*join]   # join "join" strings to path
                     [,**abs=1] # 1:  create absolute path
                                  0:  use relative path
                     )
    RETURNS: Object with obj.all
                            .dirname
                            .basename
                            .splitname
                            .splitext
    ADDITIONAL: e.g. 'C:/user/something\\new/'
                           .all    'C:\\user\\something\\new'
                           .dirname 'C:\\user\\something'
                           .basename 'new'
                           .splitname 'new'
                           .splitext  ''
                     'C:/user/something\\new/file.txt'
                           .all    'C:\\user\\something\\new\\file.txt'
                           .dirname 'C:\\user\\something\\new'
                           .basename 'file.txt'
                           .splitname 'file'
                           .splitext  '.txt'
    '''
    if abs.get('abs', 1) == 1:
        whole = os.path.normpath(
            os.path.abspath(os.path.join(pathname, *join)))
    else:
        whole = os.path.normpath(os.path.join(pathname, *join))
    path = os.path.dirname(whole)
    file_name = os.path.basename(whole)
    pref, ext = os.path.splitext(file_name)
    class NormPath:
        def __init__(self):
            self.all = whole
            self.dirname = path
            self.basename = file_name
            self.splitname = pref
            self.splitext = ext
    output = NormPath()
    return output
def bni_ramp_surf(sele, mapp, command='surface'):
    '''Color Surface by ramp.
    '''
    if command == 'mesh':
        surfname = 'mesh'
        surfcommand = 'mesh_color'
    else:
        surfname = 'surface'
        surfcommand = 'surface_color'
    print("# Color %s by map" % surfname)
    print("# --------------------")
    print("# BNI-Tools v.%s" % __version__)
    if cmd.get_type(mapp) == "object:map":
        rampname = check_names(mapp + '_ramp')
        cmd.ramp_new(rampname, mapp)
        print("# Map    : %s" % mapp)
    else:
        rampname = mapp
        mapp = rampname.split("_ramp")[0]
#    if len(filenamelist)==2:
#             cmd.set("surface_color",objramp,objpdb)
#             print "File   : %s"%objpdb
#           else:
#              print "Load Pdb file as object into pymol and"
#              print "set surface_color,%s,(object)"%objramp
#          else:
#             print "No Map loaded ... (extention has to be .phi or .dx)"
    print("# Ramp   : %s" % rampname)
    print("# %7s: %s" % (surfname, sele))
    print("--------------")
    print("To change ramp values use ctr+mid mouse button on ramp")
    print("or use")
    print("ramp_new %s,%s,[min,mid,max],[red,white,blue]" % (rampname, mapp))
    print("to set ramp values.")
    print("# %s of %s colored by" % (surfname, sele))
    print("# Ramp: %s" % rampname)
    cmd.set(surfcommand, rampname, sele)
class Bni_func_pymol:
    '''Class for passing BNI-functions in plug-in for PyMOL sidebar calls and actions.
    '''
    def __init__(self):
        self.new_views = new_views 
        self.select_hp = select_hp 
        self.get_surface_handler = get_surface_handler
        self.fast_sym = fast_sym 
        self.set_surface = set_surface 
        self.bni_ramp_surf = bni_ramp_surf 
        self.del_enabled = del_enabled 
        self.inv_enabled = inv_enabled 
        self.combine_sele = combine_sele 
        self.get_sequence = get_sequence 
        self.get_pdb_lines = get_pdb_lines 
    def test(self):
        print("Bni sidebar actions")
stored.bni_func = Bni_func_pymol()
#
# GUI TWEAK - CHANGES TO SIDEBAR MENU for bni-tools
#
from pymol import menu
dummy = lambda *args,**kwargs: []
# MENU DECORATOR
def menuappend(f):
    '''Decorator for overloading menu functions by appending to them'''
    orig = getattr(menu, f.__name__, dummy)
    wrapped = lambda *args,**kwargs: orig(*args,**kwargs) + f(*args,**kwargs)
    setattr(menu, f.__name__, wrapped)
    return wrapped
# BNIADD
@menuappend
def mol_show(self_cmd, sele):
    return [
        #BNI-Tools Transparency and Surface
        [ 0, ''           , ''                               ],
        [ 1, 'transparency'    , bni_transparency_all(self_cmd,sele) ],
        [ 1, 'surface flag'    , bni_surface(self_cmd,sele) ]]
def bni_transparency_all(self_cmd,sele):
    result=[[2,'transp:','']]
    transitems_obj=["surface","sphere","cgo","nonbonded","space","stick","ribbon","cartoon","cartoon_ring","ellipsoid"]
    transitems=transitems_obj
    if sele in [a.lower() for a in self_cmd.get_names('objects')]:
        transitems=transitems_obj
    for transitem in transitems:
        if transitem=="space":
             result.append([0,'',''])
        else:
             result.append([1,transitem,bni_transparency(self_cmd,sele,object=transitem)])
    return result
def bni_transparency(self_cmd,sele,object="cgo",label="%:"):
    result= [[2,label,'']]
    if object=="surface": 
       object=""
    else:
       object=object+"_"
    transproc=[0,10,20,30,40,50,60,70,80,90]
    for proc in transproc:
        result.append([1,str(proc)+'%','cmd.set("%stransparency",%s,"'%(object,(proc/100.0))+sele+'")'])
    return result
#BNI-Transparency END
#BNI-Surface BEGIN
def bni_surface(self_cmd,sele):
    result= [[2,"alter surface",''],
             [1,"back to standard",'stored.bni_func.get_surface_handler("'+sele+'","Standard")'],
             [0,''                ,''                                                          ],
             ]
    surf_settings=["Methionine","Exclude","Hide","separator","Include","Show"]
    surf_settings_name=['include seleno-methionine',"exclude for atom(s)","hide for atom(s)","separator","include for atom(s)","show for atom(s)"]
    for setting_name,setting in zip(surf_settings_name,surf_settings):
        if setting_name=="separator":
             result.append([ 0, ''           , ''                               ])
             continue
        result.append([1,str(setting_name).lower(),'stored.bni_func.get_surface_handler("'+sele+'","'+setting+'")'])
    return result
@menuappend
def cgo_show(self_cmd, sele):
    return [
            #BNI-Tools transparency add to cgo
            [ 1, 'transp'    , bni_transparency(self_cmd,sele,object="cgo") ]
             ]   
#BNI Tools surf colors,lable colors
if not hasattr(menu,"by_rep"):
    @menuappend
    def by_surf(self_cmd, sele,command='cmd.set("surface_color",'):
        return [[2,'Color',''                                               ],
               [ 1, 'by atom','%s-1,selection="'%command+sele+'")'                    ],
               [ 0, ''                                , ''                  ],
               [ 1, 'by map',bni_surf_map(self_cmd,sele)                    ]
               ] + all_surf_colors(self_cmd,sele,command)
    @menuappend
    def by_mesh(self_cmd, sele,command='cmd.set("mesh_color",'):
        return [[2,'Color',''                                               ],
               [ 1, 'by atom','%s-1,selection="'%command+sele+'")'                    ],
               [ 0, ''                                , ''                  ],
               [ 1, 'by map',bni_surf_map(self_cmd,sele,command='mesh')     ]
               ] + all_surf_colors(self_cmd,sele,command)
    @menuappend
    def by_label(self_cmd, sele,command='cmd.set("label_color",'):
        return [[2,'Color',''                                               ],
               [ 1, 'by atom','%s"default",selection="'%command+sele+'")'    ],
               [ 0, ''                                , ''                 ]
               ] + all_surf_colors(self_cmd,sele,command)
    @menuappend
    def by_stick(self_cmd, sele,command='cmd.set_bond("stick_color",',selekey="selection1"):
        self_cmd.set("stick_ball_color","atom")
        return [[2,'Color',''                                               ],
               [ 1, 'by atom','%s"'%'cmd.unset_bond("stick_color",'+sele+'")'    ],
               [ 0, ''                                , ''                 ]
               ] + all_surf_colors(self_cmd,sele,command,selekey)
else:
    @menuappend
    def by_rep_sub(self_cmd, rep, setting, sele):
        r = [[ 0, ''                                , ''                 ],
             [ 1, 'my colors',owncolors_surf(self_cmd,sele,command='cmd.color(')]]
        return r
@menuappend
def by_helix(self_cmd, sele,command='cmd.color('):
    sele="(ss H and ("+sele+"))"
    return [[2,'Color',''                                               ],
           ] + all_surf_colors(self_cmd,sele,command)
@menuappend
def by_sheet(self_cmd, sele,command='cmd.color('):
    sele="(ss S and ("+sele+"))"
    return [[2,'Color',''                                               ],
           ] + all_surf_colors(self_cmd,sele,command)
@menuappend
def by_loop(self_cmd, sele,command='cmd.color('):
    sele="((not (ss S+H)) and ("+sele+"))"
    return [[2,'Color',''                                               ],
           ] + all_surf_colors(self_cmd,sele,command)
def bni_surf_map(self_cmd,sele,command='surface'):
    result=[[2,'Map',''                                               ]]
    list = self_cmd.get_names("public_objects")[0:25] # keep this practical
    list = [x for x in list if (self_cmd.get_type(x)=="object:map" or self_cmd.get_type(x)=="object:")]
    for mapp in list:
        if mapp != sele:
            result.append([1,mapp,'stored.bni_func.bni_ramp_surf("'+sele+'","'+mapp+'","'+command+'")'])
    if len(result) == 1:
        result.append([1,"no map",''])
    return result
def all_surf_colors(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
        return [
        [ 1, '\\900reds'        ,reds_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\090greens'      ,greens_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\009blues'       ,blues_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\990yellows'      ,yellows_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\909magentas'    , magentas_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\099cyans'        , cyans_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\950oranges'        , oranges_surf(self_cmd, sele, command,selekey) ],   
        [ 1, '\\978tints'        ,tints_surf(self_cmd, sele, command,selekey) ],
        [ 1, '\\666grays'        ,grays_surf(self_cmd, sele, command,selekey) ],
        [ 0, '', ''],
        [ 1, '\\999my colors'        ,owncolors_surf(self_cmd, sele, command,selekey) ],
            ]
def reds_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Reds'     ,''                               ],
        [1,'\\900red','%s4,%s="'%(command,selekey)+sele+'")'],
        [1,'\\922tv_red','%s32,%s="'%(command,selekey)+sele+'")'],
        [1,'\\634raspberry','%s5268,%s="'%(command,selekey)+sele+'")'],
        [1,'\\755darksalmon','%s5280,%s="'%(command,selekey)+sele+'")'],
        [1,'\\955salmon','%s9,%s="'%(command,selekey)+sele+'")'],
        [1,'\\944deepsalmon','%s5258,%s="'%(command,selekey)+sele+'")'],
        [1,'\\824warmpink','%s5279,%s="'%(command,selekey)+sele+'")'],
        [1,'\\611firebrick','%s49,%s="'%(command,selekey)+sele+'")'],
        [1,'\\522ruby','%s21,%s="'%(command,selekey)+sele+'")'],
        [1,'\\521chocolate','%s50,%s="'%(command,selekey)+sele+'")'],
        [1,'\\632brown','%s51,%s="'%(command,selekey)+sele+'")'],
        ]
def greens_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Greens'     ,''                               ],
        [1,'\\090green','%s3,%s="'%(command,selekey)+sele+'")'],
        [1,'\\292tv_green','%s33,%s="'%(command,selekey)+sele+'")'],
        [1,'\\490chartreuse','%s14,%s="'%(command,selekey)+sele+'")'],
        [1,'\\570splitpea','%s5267,%s="'%(command,selekey)+sele+'")'],
        [1,'\\564smudge','%s5270,%s="'%(command,selekey)+sele+'")'],
        [1,'\\686palegreen','%s5259,%s="'%(command,selekey)+sele+'")'],
        [1,'\\094limegreen','%s15,%s="'%(command,selekey)+sele+'")'],
        [1,'\\494lime','%s10,%s="'%(command,selekey)+sele+'")'],
        [1,'\\792limon','%s5276,%s="'%(command,selekey)+sele+'")'],      
        [1,'\\252forest','%s22,%s="'%(command,selekey)+sele+'")'],
        ]
def blues_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Blues'     ,''                               ],
        [1,'\\009blue','%s2,%s="'%(command,selekey)+sele+'")'],
        [1,'\\339tv_blue','%s34,%s="'%(command,selekey)+sele+'")'],
        [1,'\\049marine','%s17,%s="'%(command,selekey)+sele+'")'],
        [1,'\\449slate','%s11,%s="'%(command,selekey)+sele+'")'],
        [1,'\\779lightblue','%s5263,%s="'%(command,selekey)+sele+'")'],
        [1,'\\247skyblue','%s5277,%s="'%(command,selekey)+sele+'")'],
        [1,'\\409purpleblue','%s16,%s="'%(command,selekey)+sele+'")'],
        [1,'\\226deepblue','%s23,%s="'%(command,selekey)+sele+'")'],
        [1,'\\115density','%s4155,%s="'%(command,selekey)+sele+'")'],
        ]
def yellows_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Yellows'     ,''                               ],
        [1,'\\990yellow','%s6,%s="'%(command,selekey)+sele+'")'],
        [1,'\\992tv_yellow','%s35,%s="'%(command,selekey)+sele+'")'],
        [1,'\\994paleyellow','%s5256,%s="'%(command,selekey)+sele+'")'],
        [1,'\\983yelloworange','%s36,%s="'%(command,selekey)+sele+'")'],            
        [1,'\\792limon','%s5276,%s="'%(command,selekey)+sele+'")'],
        [1,'\\976wheat','%s52,%s="'%(command,selekey)+sele+'")'],
        [1,'\\653sand','%s5269,%s="'%(command,selekey)+sele+'")'],
        ]
def magentas_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Magentas'     ,''                               ],
        [1,'\\909magenta','%s8,%s="'%(command,selekey)+sele+'")'],
        [1,'\\927lightmagenta','%s154,%s="'%(command,selekey)+sele+'")'],
        [1,'\\904hotpink','%s12,%s="'%(command,selekey)+sele+'")'],
        [1,'\\968pink','%s48,%s="'%(command,selekey)+sele+'")'],
        [1,'\\978lightpink','%s5274,%s="'%(command,selekey)+sele+'")'],
        [1,'\\644dirtyviolet','%s5272,%s="'%(command,selekey)+sele+'")'],
        [1,'\\949violet','%s53,%s="'%(command,selekey)+sele+'")'],
        [1,'\\525violetpurple','%s5271,%s="'%(command,selekey)+sele+'")'],
        [1,'\\707purple','%s19,%s="'%(command,selekey)+sele+'")'],
        [1,'\\515deeppurple','%s5261,%s="'%(command,selekey)+sele+'")'],
        ]
def cyans_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Cyans'     ,''                               ],
        [1,'\\099cyan','%s5,%s="'%(command,selekey)+sele+'")'],
        [1,'\\799palecyan','%s5265,%s="'%(command,selekey)+sele+'")'],
        [1,'\\499aquamarine','%s5257,%s="'%(command,selekey)+sele+'")'],
        [1,'\\297greencyan','%s5275,%s="'%(command,selekey)+sele+'")'],
        [1,'\\077teal','%s20,%s="'%(command,selekey)+sele+'")'],
        [1,'\\155deepteal','%s5262,%s="'%(command,selekey)+sele+'")'],
        [1,'\\466lightteal','%s5266,%s="'%(command,selekey)+sele+'")'],
        ]
def oranges_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Oranges'     ,''                               ],
        [1,'\\950orange','%s13,%s="'%(command,selekey)+sele+'")'],
        [1,'\\951tv_orange','%s37,%s="'%(command,selekey)+sele+'")'],
        [1,'\\962brightorange','%s30,%s="'%(command,selekey)+sele+'")'],
        [1,'\\985lightorange','%s5264,%s="'%(command,selekey)+sele+'")'],      
        [1,'\\983yelloworange','%s36,%s="'%(command,selekey)+sele+'")'],      
        [1,'\\760olive','%s18,%s="'%(command,selekey)+sele+'")'],
        [1,'\\551deepolive','%s5260,%s="'%(command,selekey)+sele+'")'],
        ]
def tints_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Tints'     ,''                               ],
        [1,'\\976wheat','%s52,%s="'%(command,selekey)+sele+'")'],
        [1,'\\686palegreen','%s5259,%s="'%(command,selekey)+sele+'")'],
        [1,'\\779lightblue','%s5263,%s="'%(command,selekey)+sele+'")'],      
        [1,'\\994paleyellow','%s5256,%s="'%(command,selekey)+sele+'")'],
        [1,'\\978lightpink','%s5274,%s="'%(command,selekey)+sele+'")'],
        [1,'\\799palecyan','%s5265,%s="'%(command,selekey)+sele+'")'],
        [1,'\\985lightorange','%s5264,%s="'%(command,selekey)+sele+'")'],            
        [1,'\\889bluewhite','%s5278,%s="'%(command,selekey)+sele+'")'],
        ]
def grays_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    return [
        [ 2, 'Grays'     ,''                               ],
        [ 1, '\\999white ', '%s"white",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\999gray90 ', '%s"grey90",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\888gray80 ', '%s"grey80",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\777gray70 ', '%s"grey70",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\666gray60 ', '%s"grey60",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\555gray50 ', '%s"grey50",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\444gray40 ', '%s"grey40",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\333gray30 ', '%s"grey30",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\222gray20 ', '%s"grey20",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\222gray10 ', '%s"grey10",%s="'%(command,selekey)+sele+'")'  ],
        [ 1, '\\222black ', '%s"black",%s="'%(command,selekey)+sele+'")'  ],
        ]
def owncolors_surf(self_cmd, sele,command='cmd.set("surface_color",',selekey="selection"):
    owncolors=[(colname,colnumber) for (colname,colnumber) in self_cmd.get_color_indices() if int(colnumber)>=(5387)]
    result=[[ 2, 'My Colors'     ,''                               ]]
    for (colname,colnumber) in owncolors:
        result.append([ 1, '\\999%s '%colname, '%s%s,%s="'%(command,colnumber,selekey)+sele+'")'  ])
    return result
@menuappend
def by_ss(self_cmd, sele):
    return [
    #BNI-Tools color
    [ 0, ''                                    , ''                      ],
    [ 1, '\\900H\\950e\\990l\\090i\\099x'      , by_helix(self_cmd, sele)],
    [ 1, '\\900S\\950c\\990h\\090e\\099e\\059t', by_sheet(self_cmd, sele)],
    [ 1, '\\900L\\950o\\990o\\090p'            , by_loop(self_cmd, sele) ]
    ]
#BNI-Surface END
@menuappend
def spectrum(self_cmd, sele):
    r = [
        #BNI-Tools color by occupancy
        [ 0, ''                                , ''                 ],
        [ 1, 'occupancy'   , 'cmd.spectrum("q",selection=("'+sele+'"),quiet=0)'         ],
        [ 1, 'occupancy(*/ca)'   , 'cmd.spectrum("q",selection="(('+sele+')&*/ca)",quiet=0)'         ],
        ]
    return r
# check for older menu versions of PyMOL
# check if by_rep is present
if hasattr(menu,"menucontext") and hasattr(menu,"by_rep"):
    @menuappend
    def mol_color(self_cmd, sele):
        with menu.menucontext(self_cmd, sele):
          return (
            [[ 0, ''                                , ''                 ],
             [ 1, 'my colors',owncolors_surf(self_cmd,sele,command='cmd.color(')],])
elif hasattr(menu,"menucontext"):
    @menuappend
    def mol_color(self_cmd, sele):
        with menu.menucontext(self_cmd, sele):
          return (
            [[ 0, ''                                , ''                 ],
             [ 1, 'my colors',owncolors_surf(self_cmd,sele,command='cmd.color(')],
             [ 0, ''                                , ''                 ],
             [ 1, 'surface'  ,by_surf(self_cmd, sele) ],
             [ 1, 'mesh'     ,by_mesh(self_cmd, sele) ],
             [ 1, 'label'    ,by_label(self_cmd,sele) ],
             [ 1, 'stick'    ,by_stick(self_cmd,sele) ],
             [ 0, ''                                , ''                 ],
             ] )
else:
    @menuappend
    def mol_color(self_cmd, sele):
          return (
            [[ 0, ''                                , ''                 ],
             [ 1, 'my colors',owncolors_surf(self_cmd,sele,command='cmd.color(')],
             [ 0, ''                                , ''                 ],
             [ 1, 'surface'  ,by_surf(self_cmd, sele) ],
             [ 1, 'mesh'     ,by_mesh(self_cmd, sele) ],
             [ 1, 'label'    ,by_label(self_cmd,sele) ],
             [ 1, 'stick'    ,by_stick(self_cmd,sele) ],
             [ 0, ''                                , ''                 ],
             ] )
def bni_preset_hp(self_cmd,sele):
    return [[ 2, 'HP values:', ''],
              [ 1, 'KandD'   , 'stored.bni_func.select_hp("KandD","'+sele+'")'          ],
              [ 1, 'Rose'    , 'stored.bni_func.select_hp("Rose","'+sele+'")'          ],
              [ 1, 'GES'     , 'stored.bni_func.select_hp("GES","'+sele+'")'           ],
              ]
def bni_preset_sym(self_cmd,sele):
    result=[[ 2, 'A', '']]
    distances=[2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5]
    for dist in distances:
       result.append([ 1, str(dist)   , 'stored.bni_func.fast_sym(cutoff="'+str(dist)+'",selection="'+sele+'")'])
    return result
def bni_preset_symops(self_cmd,sele):
    return [[ 2, 'show symmetry patch within',''],
            [ 1, 'set distance',bni_preset_sym(self_cmd,sele)],
            [ 1, 'set surface patch color',bni_preset_sym_color(self_cmd,sele)]]
def bni_preset_sym_color(self_cmd,sele):
    return [[ 2, 'surface patch color',''],
            [ 1, 'to symmetry representation','stored.bni_func.set_surface(mode="symmetry")'],
            [ 1, 'to standard representation','stored.bni_func.set_surface(mode="standard")']]
#BNI Presets end
@menuappend
def presets(self_cmd, sele):
    return [
              #BNI-Tools 0.3 BEGIN
              [ 0, ''               ,''                             ],
              [ 1, 'track main chain'   ,'stored.bni_func.new_views("track_mainchain","'+sele+'")'],
              [ 1, 'symmetry surface'   ,bni_preset_symops(self_cmd,sele)],
              [ 1, 'hydrophobic residues'   , bni_preset_hp(self_cmd,sele)         ],
              [ 1, 'surface inspection'   , 'stored.bni_func.new_views("surf_rel_res","'+sele+'")'],
              ]
@menuappend
def sele_action(self_cmd, sele):
    return [
              #BNI-Tools sequence begin
              [ 0, ''          ,''                                              ],
              [ 1, 'sequences', bni_sequence(self_cmd,sele)],
              [ 0, ''          ,''                                              ],
              #BNI Tools sequence end
              ]
@menuappend
def group_action(self_cmd, sele):
    return [
            #BNI-Tools sequence begin
            [ 0, ''          ,''                                              ],
            [ 1, 'sequences', bni_sequence(self_cmd,sele)],
            [ 0, ''          ,''                                              ],
            #BNI-Tools sequence end
            ]
@menuappend
def mol_action(self_cmd, sele):
    return [
            #BNI-Tools sequence begin
            [ 0, ''          ,''                                              ],
            [ 1, 'sequences', bni_sequence(self_cmd,sele)],
            [ 0, ''          ,''                                              ],
            #BNI-Tools sequence end
              ]
@menuappend
def all_action(self_cmd, sele):
    return [
              #BNI-Tools BEGIN
              [ 0, ''          ,''                                              ],
              [ 1, 'delete enabled', bni_selection_del(self_cmd,sele)],
              [ 1, 'invert enabled/disabled', bni_selection_switch(self_cmd,sele)],
              [ 1, 'combine selections',bni_selection_combine(self_cmd,sele)],
              [ 0, ''          ,''                                              ],
              ]
#BNI-Selection BEGIN 
def bni_sequence(self_cmd,sele):
    return [[2,'Format',''],
            [1,'fasta','stored.bni_func.get_sequence(selection="'+sele+'",format="fasta")'],
            [1,'pir','stored.bni_func.get_sequence(selection="'+sele+'",format="pir")'],
            [1,'modeller','stored.bni_func.get_sequence(selection="'+sele+'",format="modeller")'],
            [1,'list',bni_sequence_list(self_cmd,sele)],
            [1,'pdb','stored.bni_func.get_pdb_lines(selection="'+sele+'",state="0",enabled="1")']]
def bni_sequence_list(self_cmd,sele):
   return [[ 2, 'List',''],
           [ 1, 'residues','stored.bni_func.get_sequence(selection="'+sele+'",format="normal")'],
           [ 1, 'bfac per res',bni_sequence_all(self_cmd,sele,format="normalb")],
           [ 1, 'occ  per res',bni_sequence_all(self_cmd,sele,format="normalq")],
           [ 1, 'all atoms','stored.bni_func.get_sequence(selection="'+sele+'",selecting="atom",format="normal")']]
def bni_sequence_all(self_cmd,sele,format):
   return [[ 2,'include',''],
	       [ 1,'all atoms','stored.bni_func.get_sequence(selection="'+sele+'",format="%s_all")'%format],
		   [ 1,'mainchain','stored.bni_func.get_sequence(selection="'+sele+'",format="%s_main")'%format],
	       [ 1,'sidechain','stored.bni_func.get_sequence(selection="'+sele+'",format="%s_side")'%format]]
def bni_selection_switch(self_cmd,sele):
    return [[2,'Invert',''],
            [1,'all','stored.bni_func.inv_enabled("all")'],
            [1,'objects','stored.bni_func.inv_enabled("obj")'],
            [1,'selecions','stored.bni_func.inv_enabled("selections")']]
def bni_selection_combine(self_cmd,sele):
    return [[2,'Combine to (sele)',''],
            [1,'all','stored.bni_func.combine_sele(check="all")'],
            [1,'enabled','stored.bni_func.combine_sele(check="enabled")'],
            [1,'disabled','stored.bni_func.combine_sele(check="disabled")']]
def bni_selection_del(self_cmd,sele):
    return [[2,'Delete Enabled',''],
            [1,'all','stored.bni_func.del_enabled("all")'],
            [1,'objects','stored.bni_func.del_enabled("obj")'],
            [1,'selecions','stored.bni_func.del_enabled("selections")']]
#BNI-Selection END
# BNIADDEND
