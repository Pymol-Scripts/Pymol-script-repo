from pymol.controlling import ring_dict,mode_name_dict,mode_dict
 
# redefine the three_button_viewing mode
mode_name_dict['three_button_viewing'] = 'My 3-But View'
mode_dict['three_button_viewing'] =  [ ('l','none','rota'),
                      ('m','none','move'),
                      ('r','none','movz'),
                      ('l','shft','+Box'),
                      ('m','shft','-Box'),
                      ('r','shft','clip'),                 
                      ('l','ctrl','+/-'),
                      ('m','ctrl','pkat'),
                      ('r','ctrl','pk1'),                 
                      ('l','ctsh','Sele'),
                      ('m','ctsh','orig'),
                      ('r','ctsh','menu'),
                      ('w','none','slab'),
                      ('w','shft','movs'),
                      ('w','ctrl','mvsz'),
                      ('w','ctsh','movz'),
                      ('double_left','none','menu'),
                      ('double_middle','none','none'),
                      ('double_right','none', 'pk1'),
                      ('single_left','none','+/-'),
                      ('single_middle','none','cent'),
                      ('single_right','none', 'pkat'),
                      ]
