reinitialize
import elbow_angle

bg_color white

fetch 3ghe, async=0
dss
as cartoon
util.cbc

hide everything, all and not polymer

set_view (\
     0.155146211,    0.806551158,    0.570442557,\
     0.960532069,    0.011798559,   -0.277920276,\
    -0.230887875,    0.591045439,   -0.772885919,\
     0.000000000,    0.000000000, -276.763153076,\
    18.085998535,   -4.697593689,   16.178558350,\
   226.579132080,  326.947174072,  -20.000000000 )

# use defaults: light=L, heavy=H, limit_l=107, limit_h=113
elbow_angle 3ghe, draw=1