reinitialize 
fetch 1HXB, async=0
#load 1hxb.pdb, 1HXB
create HIV-1_PROTEINASE, 1HXB and polymer and not resn ROC
create ROC_A, 1HXB and resn ROC and alt a
create ROC_B, 1HXB and resn ROC and alt b
delete 1HXB
 
hide everything, all
#h_add cdk2
#h_add EFP
show_as cartoon, HIV-1_PROTEINASE
show_as sticks, ROC_A or ROC_B
util.cbay ROC_A
util.cbap ROC_B
 
select flexible_A, byres HIV-1_PROTEINASE within 3.5 of ROC_A
select flexible_B, byres HIV-1_PROTEINASE within 3.5 of ROC_B
show sticks, flexible_A or flexible_B
util.cbag flexible_A or flexible_B 
disable flexible_A or flexible_B
