reinitialize
import elbow_angle

bg_color white

# example structures from Stanfield, et al. JMB 2006
# doi:10.1016/j.jmb.2006.01.023
fetch 1bbd, async=0
fetch 7fab, async=0
fetch 1dba, async=0
fetch 1plg, async=0
fetch 1nl0, async=0

dss
as cartoon
set cartoon_transparency, 0.7

remove all and not chain L+H
util.mass_align("1bbd and ((chain L and resi 1-114) or (chain H and resi 1-118))")

# adopt a similar view to Figure 1 in Stanfield, et al.
set_view (\
    -0.953261435,   -0.226005003,    0.200535893,\
    -0.230026290,    0.112494141,   -0.966659248,\
     0.195909262,   -0.967606425,   -0.159222543,\
     0.000000000,    0.000000000, -230.122619629,\
    62.279075623,   48.879341125,  138.177505493,\
   181.430419922,  278.814819336,  -20.000000000 )

# 1bbd
# Stanfield:      127 deg
# elbow_angle.py: 125 deg
elbow_angle 1bbd, limit_l=114, limit_h=118, draw=1

# 7fab
# Stanfield:      132 deg
# elbow_angle.py: 126 deg
elbow_angle 7fab, limit_l=104, limit_h=117, draw=1

# 1dba
# Stanfield:      183 deg
# elbow_angle.py: 176 deg
elbow_angle 1dba, draw=1

# 1plg
# Stanfield:      190 deg
# elbow_angle.py: 189 deg
elbow_angle 1plg, limit_l=112, limit_h=117, draw=1

# 1nl0
# Stanfield:      220 deg
# elbow_angle.py: 221 deg
elbow_angle 1nl0, draw=1