def getChargeMass(atoms):
    elist = []
    totmass = 0.0
    
    for a in atoms:
        if a.element=='H':
            e=1.0
            totmass+=1.00794
        elif a.element=='C':
            e=6.0
            totmass+=12.0107
        elif a.element=='AU':
            e=79.0
            totmass+=196.96655
        elif a.element=='N':
            e=7.0
            totmass+=14.00674
        elif a.element=='O':
            e=8.0
            totmass+=15.9994
        elif a.element=='P':
            e=15.0
            totmass+=30.973761
        elif a.element=='S':
            e=16.0
            totmass+=32.066
        elif a.element=='W':
            e=18.0
            totmass+=1.00794*2+15.9994 # ficticious water 'atom'
        else:
            print "skipping unknown atom %s"%a.element
            e = -1
        elist.append(e)

    return elist, totmass

