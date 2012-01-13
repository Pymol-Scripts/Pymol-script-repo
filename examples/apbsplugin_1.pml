reinitialize 
fetch 3IG7, async=0
#load 3ig7.pdb, 3IG7
create cdk2, 3IG7 and polymer
create EFP, 3IG7 and organic and not resn ACE
delete 3IG7
 
hide everything, all
#h_add cdk2
#h_add EFP
show_as cartoon, cdk2
show_as sticks, EFP
util.cbay EFP
 
select flexible, byres cdk2 within 3.5 of EFP
show sticks, flexible
util.cbag flexible
disable flexible
zoom cdk2
