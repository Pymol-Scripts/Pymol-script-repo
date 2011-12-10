##
##Authors : Michel F Sanner , Sowjanya Karnati
##
##
import sys,getopt,os
import types,string

######################################################################
#   COMMAND LINE OPTIONS
######################################################################
def Usage():
    "Print helpful accurate statement to Usage statement"
    print "Usage:getlatest"
    print "Description of Command:"
    print "-----------------------"
    print "This command cvs co and installs package or pacakges "
    print "Optional parameters:"
    print "--------------------"
    print "-p   :   packages seperated by commas or a pacakge"
    print "         eg: ./getlatest -psymserv,gle,opengltk "
    print "         or ./getlatest -psymserv" 
    print "         or  ./getlatest -pall ,it updates all the packages "
    print "         or ./getlatest -p platdep ,Updates only platform dependent packages ie, eg: opengltk,gle,bhtree etc"
    print  "        or ./getlatest -p platindep ,Updates only platform independent packages ie, pure python packages eg: Pmv,Vision,AutoDockTools etc"
    print "-t   :   tag to use with cvs"
    print "         [consistent with all other packages] " 
    print "-h   :    help"


try:
    exec("opt_list,args=getopt.getopt(sys.argv[1:],'p:t:h')")
except getopt.GetoptError,msg:
    print "getlatest:'%s'"%msg
    Usage()
    sys.exit(2)

#initialize required parameters
package = None

cvsTag = None

for o,a in opt_list:
    if o in ('-p','--p'):
        package = a
    if o in ('-t', '--t'):
        cvsTag = a
    if o in ('-h','--h'):
        Usage()
        sys.exit()

###########MGLPACKS ##########################



DPACKAGES =['bhtree','gle','geomutils','mslib','opengltk','sff','stride','UTpackages','QSlimLib','cAutoDock','OpenCSG']

INDPACKAGES =['AutoDockTools','Vision','symserv','DejaVu','NetworkEditor','ViewerFramework','Pmv','PyAutoDock','mglutil','Volume','MolKit','PyBabel']

DISTPACKAGES = DPACKAGES  + INDPACKAGES


####caluculating todays date
datestr = os.system("date '+%m/%d/%y' >datefile.txt")
fptr=open('datefile.txt','r')
datelines=fptr.readlines()
todays_date = str(datelines[0]).replace('/','-')[:-1]
os.system('rm -rf datefile.txt')

if not package:
    print "packupdate requires package name to be updated"
    sys.exit()

#find $PATH
#os.popen("echo $PATH >myfile.txt")    
#fptr=open('myfile.txt')
#allines= fptr.readlines()
#x = allines[0].replace(':',',')
#pathstring =x.split(',')
#fptr.close()
#os.system('rm -rf myfile.txt')

pathstring=string.split(os.environ["PATH"], ":") 
############FIND SWIG#################


def findswig(pathstring):
    list = []
    from distutils.command import build_ext
    from distutils.dist import Distribution
    d  =Distribution()
    a = build_ext.build_ext(d)
    yourswig = a.find_swig()
    
    for a in pathstring:
        if  os.path.exists('%s/%s' %(a,yourswig)):
            list.append('%s/%s' %(a,yourswig))  
             
    if list == []:
        print "YOU DON'T HAVE SWIG"
        sys.exit()
        
############FIND CC##################

def findcc(pathstring):
    list = []
    from distutils import ccompiler    
    yourcompiler = ccompiler.get_default_compiler()
    if sys.platform == 'sunos5"' or 'irix6':
        for a in pathstring:
            if  os.path.exists('%s/cc' %a) or os.path.exists('%s/CC' %a):
                list.append('%s' %a)    
        if list == []:
                print "YOU DO NOT HAVE A C-COMPILER"
                sys.exit()
    if (sys.platform == 'cygwin') or (os.name == 'nt'):
     	for a in pathstring:
            if  os.path.exists('%s/cc' %a):
                list.append('%s' %a)    
        if list == []:
                print "YOU DO NOT HAVE A C-COMPILER"
                sys.exit()		

    if (sys.platform == 'linux2') or (sys.platform == 'darwin'):
        for a in pathstring:
            if  os.path.exists('%s/gcc' %a) or os.path.exists('%s/g++' %a):
                list.append('%s' %a)    
        	
	    if list == []:
                print "YOU DO NOT HAVE A C-COMPILER"
                sys.exit()

#############FIND CVS#################
def findcvs(pathstring):
    list =[]
    for a in pathstring:
        if  os.path.exists('%s/cvs' %a):
            list.append('%s/cvs' %a)
    if list == []:
        print "YOU DO NOT HAVE CVS"
        sys.exit()

                
    
        
findswig(pathstring)
findcc(pathstring)
findcvs(pathstring)


   
uninstlist =[]
curdir = os.getcwd()

curdir=  os.environ["MGL_ROOT"]

archosvdir = os.environ["MGL_ARCHOSV"]
py_version =  (string.split(sys.version))[0][0:3]

dpdir=os.path.join(curdir, archosvdir, "lib", "python"+py_version, "site-packages")
indpdir=os.path.join(curdir,"share", "lib", "python"+py_version, "site-packages")
bindir=os.path.join(curdir ,archosvdir, "bin")

if os.path.exists(os.path.join(dpdir, "MGLPACKS")):
    os.system('rm -rf %s' % os.path.join(dpdir, "MGLPACKS"))

PACKLIST=[]
uninstlist=[]
if package=='all':
    os.chdir(dpdir)    
    DPLIST=os.listdir('./')
    os.chdir(indpdir)    
    INDPLIST=os.listdir('./')
    PACKLIST=DPLIST+INDPLIST
    PACK_DIST=[]
    for p in PACKLIST:
        if p in DISTPACKAGES:
            PACK_DIST.append(p)
    	 

if package=='platdep':
    os.chdir(dpdir)    
    PACKLIST=os.listdir('./')    
    PACK_DIST=[]
    for p in PACKLIST:
        if p in DISTPACKAGES:
            PACK_DIST.append(p)

if package=='platindep':
    os.chdir(indpdir)    
    PACKLIST=os.listdir('./')    
    PACK_DIST=[]
    for p in PACKLIST:
        if p in DISTPACKAGES:
            PACK_DIST.append(p)

if not package in ['all','platdep','platindep']:
    PACK_DIST=[]
    PACKLIST = string.split(package,',')
    for p in PACKLIST:
       PACK_DIST.append(p)
       
print "PACK_DIST:", PACK_DIST

if PACK_DIST!=[]:
    for p in PACK_DIST:
        if p in DPACKAGES:
            if not os.path.exists('%s/MGLPACKS' %(dpdir)):
                os.mkdir('%s/MGLPACKS' %dpdir)
            if not os.path.exists('%s/MGLPACKS%s' %(dpdir,todays_date)):
                os.mkdir('%s/MGLPACKS%s' %(dpdir,todays_date))
        if p in INDPACKAGES:
            if not os.path.exists('%s/MGLPACKS%s' %(indpdir,todays_date)):
                os.mkdir('%s/MGLPACKS%s' %(indpdir,todays_date))

        if p in PACKLIST and p in DPACKAGES:
            os.chdir('%s' %dpdir)
            if not os.path.exists('./MGLPACKS%s/%s' %(todays_date,p)):
                os.system('mv %s MGLPACKS%s/%s' %(p,todays_date,p))
            else:
                os.system('rm -rf %s' %p)
        if p in PACKLIST and p in INDPACKAGES:
            os.chdir(indpdir)
            if not os.path.exists('./MGLPACKS%s/%s' %(todays_date,p)): 
                os.system('mv %s MGLPACKS%s/%s' %(p,todays_date,p))
            else:
                os.system('rm -rf %s' %p)       

  

#Connecting to CVS
#CVS login
print "PRESS ENTER WHEN ASKED FOR PASSWORD"
os.system('/usr/bin/cvs -d:pserver:anonymous@moses.scripps.edu:/export/cvs login')


print PACK_DIST

if PACK_DIST==[]:
	print "No Packages are found in the current installation to update.Use getlatest -ppackname to update packages you want"
	sys.exit()

#CVS check out 
PYTHON = sys.executable
for pack in  PACK_DIST:
    if pack in DPACKAGES:
        pack=pack+"DIST"
        os.chdir(dpdir)
        if os.path.exists('./MGLPACKS'):
            os.chdir('./MGLPACKS')
        if cvsTag:
            st =os.system('/usr/bin/cvs -z3 -d:pserver:anonymous@moses.scripps.edu:/export/cvs co -r %s %s' % (cvsTag, pack))
        else:
            st =os.system('/usr/bin/cvs -z3 -d:pserver:anonymous@moses.scripps.edu:/export/cvs co -A %s' %pack)    
        if st!=0:
            print "cannot cvs co package %s" % pack
            sys.exit()
    
        os.chdir('./%s' %pack)
        #cmd = 'python2.5 setup.py install --install-platlib=%s --install-purelib=%s --install-scripts=%s'%(dpdir,indpdir,bindir)      
        #cmd = 'python setup.py install'
        cmd = "%s setup.py install --install-platlib=%s" % (PYTHON, dpdir)
        st = os.system(cmd)
        if st!=0:
            print "%s not installed" %pack
            uninstlist.append(pack)
            continue	
        else:
            print "Updating %s package success" %pack
    
    else:
        if pack in INDPACKAGES:
            os.chdir(indpdir)
            if cvsTag:
                st =os.system('/usr/bin/cvs -z3 -d:pserver:anonymous@moses.scripps.edu:/export/cvs co -r %s %s' % (cvsTag, pack))
            else:
                st =os.system('/usr/bin/cvs -z3 -d:pserver:anonymous@moses.scripps.edu:/export/cvs co -A %s' %pack)    
            if st!=0:
                print "cannot cvs co pack %s" % pack
                sys.exit()
     

if uninstlist!=[]:
	print "%s packages are not installed" %uninstlist

else:
    if os.path.exists('%s/MGLPACKS' %dpdir):
        os.chdir(dpdir)
        os.system('rm -rf MGLPACKS')
    print "Installation Success"



