#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/__init__.py,v 1.38.2.4 2009/05/15 17:11:30 rhuey Exp $
#
# $Id: __init__.py,v 1.38.2.4 2009/05/15 17:11:30 rhuey Exp $

import os, types
import sys
import getopt
import time
import socket
from string import split

from Support.version import __version__
from mglutil import __revision__

# create hostDict with hostMacros accessible by anyone
from AutoDockTools.adthosts import hostMacros
from AutoDockTools.autodockHosts import AutoDockHosts

hostDict = AutoDockHosts(hostMacros)

h= socket.gethostname()
hostDict[h]=hostDict['localhost']
hostDict[h]['host']=h
del hostDict['localhost']

# try to extend that dictionary with user specific hostMacros

# first try to find a adthost file in current directory
if os.path.isfile('./adthosts.py'):
    execfile('./adthosts.py')
    if globals().has_key('hostMacros'):
        hostDict.update(hostMacros)
elif sys.platform!='win32':
    # try to find the user's home directory
    import posix
    if 'HOME' in posix.environ.keys():
        try:
            execfile(os.path.join(posix.environ['HOME'],'adthosts.py'))
            if globals().has_key('hostMacros'):
                hostDict.update(hostMacros)
        except:
            pass
   

def setADTmode(modeStr, mv):
    import Tkinter
    from AutoDockTools.autotorsCommands import AdtSetMode
    mv.addCommand(AdtSetMode(),'ADTSetMode')
    if not hasattr(mv.GUI, 'adtBar'):
        mv.GUI.adtBar = mv.GUI.menuBars['AutoToolsBar']
    if not hasattr(mv.GUI, 'adtFrame'):
        mv.GUI.adtFrame = vf.GUI.adtBar.menubuttons.values()[0].master
    if not hasattr(mv.GUI, 'adt4ModeLabel'):
        mv.GUI.adt4ModeLabel=Tkinter.Label(mv.GUI.adtFrame, text="ADT4.0", width=6,
                             relief='sunken', borderwidth=1, fg='DarkGreen',
                             bg = 'ivory',anchor='w' )
        mv.GUI.adt4ModeLabel.pack(side='left')
    mv.GUI.adt4ModeLabel.bind("<Double-Button-1>", mv.ADTSetMode.guiCallback)
    return mv.ADTSetMode(modeStr)
     

def setdmode(mode, mv):
    """
load display commands for mode and set them as default command for new molecule
"""
    if mode=='cpk':
        mv.browseCommands('displayCommands', commands=['displayCPK'],
              log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.displayCPK)

    elif mode=='lines':
        mv.browseCommands('bondsCommands',
              commands=('buildBondsByDistance',), log=0)
        mv.addOnAddObjectCmd(mv.buildBondsByDistance)
        mv.browseCommands('displayCommands', commands=('ribbon',), log=0)
        mv.addOnAddObjectCmd(mv.displayLines)

    elif mode=='ss':
        mv.browseCommands('secondaryStructureCommands',
              commands=('ribbon',), log=0)
        mv.addOnAddObjectCmd(mv.ribbon)

    elif mode=='sb':
        mv.browseCommands('bondsCommands',
              commands=('buildBondsByDistance',), log=0)
        mv.addOnAddObjectCmd(mv.buildBondsByDistance)
        mv.browseCommands('displayCommands',
              commands=('displaySticksAndBalls',), log=0)
        mv.addOnAddObjectCmd(mv.displaySticksAndBalls, (),
                 {'cquality':8, 'bquality':10})

    elif mode=='lic':
        mv.browseCommands('bondsCommands',
                          commands=('buildBondsByDistance',), log=0)
        mv.addOnAddObjectCmd(mv.buildBondsByDistance)
        mv.browseCommands('displayCommands',
                          commands=('displaySticksAndBalls',), log=0)
        mv.addOnAddObjectCmd(mv.displaySticksAndBalls, (),
                             {'cquality':16, 'bquality':20, 'cradius':.2, 'bRad':.2})

    elif mode=='ms':
        mv.browseCommands('msmsCommands', commands=('computeMSMS',), log=0)
        mv.browseCommands('msmsCommands', commands=('displayMSMS',), log=0)
        mv.addOnAddObjectCmd(mv.computeMSMS, (), {'density':3.0})
        mv.addOnAddObjectCmd(mv.displayMSMS)
    
    elif mode=='ca':
        mv.browseCommands('traceCommands', commands=('computeTrace',), log=0)
        mv.browseCommands('traceCommands', commands=('extrudeTrace',), log=0)
        mv.browseCommands('traceCommands', commands=('displayTrace',), log=0)
        mv.addOnAddObjectCmd(mv.computeTrace)
        mv.addOnAddObjectCmd(mv.extrudeTrace)
        mv.addOnAddObjectCmd(mv.displayTrace)

    elif mode=='bt':
        mv.browseCommands('bondsCommands',
              commands=('buildBondsByDistance',), log=0)
        mv.browseCommands('displayCommands',
              commands=('displayBackboneTrace',), log=0)
        mv.addOnAddObjectCmd(mv.displayBackboneTrace, (), 
                 {'cquality':8, 'bquality':10, 'cradius':0.25,
                  'bRad':0.33} )
        
    elif mode=='sp':
        mv.browseCommands('splineCommands', commands=('computeSpline',), log=0)
        mv.browseCommands('splineCommands', commands=('extrudeSpline',), log=0)
        mv.browseCommands('splineCommands',
              commands=('displayExtrudedSpline',), log=0)
        mv.addOnAddObjectCmd(mv.computeSpline)
        mv.addOnAddObjectCmd(mv.extrudeSpline)
        mv.addOnAddObjectCmd(mv.displayExtrudedSpline)

    elif mode=='sssb':
        mv.browseCommands('displayCommands',
              commands=('displaySSSB',), log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.displaySSSB)


def setcmode(mode, mv):
    """
load color commands for mode and set them as default command for new molecule
"""
    if mode=='ca':
        mv.browseCommands('colorCommands', commands=('colorByAtomType',),
              log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.colorByAtomType)

    elif mode=='cr':
        mv.browseCommands('colorCommands',
              commands=('colorByResidueType',), log=0)
        mv.addOnAddObjectCmd(mv.colorByResidueType)
        
    elif mode=='cc':
        mv.browseCommands('colorCommands', commands=('colorByChains',),
              log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.colorByChains)

    elif mode=='cm':
        mv.browseCommands('colorCommands', commands=('colorByMolecules',),
                   log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.colorByMolecules)

    elif mode=='cdg':
        mv.browseCommands('colorCommands', commands=('colorAtomsUsingDG',),
              log=0, package='Pmv')
        mv.addOnAddObjectCmd(mv.colorAtomsUsingDG)

    elif mode=='cs':
        mv.browseCommands('colorCommands',
              commands=('colorResiduesUsingShapely',), log=0 )
        mv.addOnAddObjectCmd(mv.colorResiduesUsingShapely)

    elif mode=='css':
        mv.browseCommands('secondaryStructureCommands',
              commands=('colorBySecondaryStructure',), log=0)
        mv.addOnAddObjectCmd(mv.colorBySecondaryStructure)


def runADT(*argv, **kw):
    """The main function for running AutoDockTools
"""
    import sys

    if type(argv) is tuple:
        if len(argv) == 0:
            argv = None
        elif len(argv) == 1:
            argv = argv[0]
            if type(argv) is not list:
                argv = [argv]
        else:
            argv = list(argv)
    if kw.has_key("ownInterpreter"):
        ownInterpreter = kw["ownInterpreter"]
    else:
        if argv is None:
            argv = ['AutoDockTools/bin/runADT.py', '-i']
            ownInterpreter = False
        elif argv[0].endswith('runADT.py') is False:
            argv.insert(0,'-i')
            argv.insert(0,'AutoDockTools/bin/runADT.py')
            ownInterpreter = False
        else:
            ownInterpreter = True

    optlist, args = getopt.getopt(argv[1:], 'haid:c:v:', [
                'help', 'again', 'overwriteLog', 'uniqueLog', 'noLog', 'die',
            'customizer=', 'interactive', 'dmode=', 'cmode=', 'noSplash', 'vision'] )
    
    help_msg = """usage: adt <options>
            -h or --help       : print this message
            -a or --again      : play back lastlog file
            --overwriteLog     : overwrite log file
            --uniqueLog        : create a log file with a unique name
            --noLog            : turn off logging
            --noSplash         : turn off Splash Screen
            --die              : do not start GUI event loop
            --customizer file  : run the user specified file
            --lib packageName  : add a libraries of commands
            -v r or --vision run  : run vision networks on the command line
            -v o or --vision once : run vision networks and exit ADT

        -d or --dmode modes : specify a display mode
                modes can be any a combination of display mode
               'cpk'  : cpk
               'lines': lines
               'ss'   : secondary structure ribbon
               'sb'   : sticks and balls
               'lic'  : licorice
               'ms'   : molecular surface
               'ca'   : C-alpha trace
               'bt'   : backbone trace
               'sp'   : CA-spline
               'sssb' : secondary structure for proteins,
                        sticks and balls for other residues with bonds
                    lines for other residues without bonds
    
        -c or --cmode modes : specify a dispaly mode
                    color scheme:
                'ca' : color by atom
                'cr' : color by residue (RASMOL scheme)
                'cc' : color by chain
                'cm' : color by molecule
                'cdg': color using David Goodsell's scheme
                'cs' : color residues using Shapely scheme
                'css': color by secondary structure element
    
              example:
              display protein as ribbon, non protein as sticks and balls
              and color by atom type
                 adt -i --dmode sssb --cmode cr myprot.pdb
                 adt -i -m sssb -c cr myprot.pdb
    
    """
    
    customizer = None
    logmode = 'overwrite'
    libraries = []
    again = 0
    interactive = 0
    die=0
    noSplash = False
    dmode = cmode = None
    dmodes = ['cpk', 'lines', 'ss', 'sb', 'lic', 'ms', 'ca', 'bt', 'sp', 'sssb' ]
    cmodes = ['ca', 'cr', 'cc', 'cm', 'cdg', 'cs', 'css']
    visionarg = None

    for opt in optlist:
        if opt[ 0] in ('-h', '--help'):
            print help_msg
            sys.exit()
        elif opt[ 0] in ('-a', '--again'): 
            again = 1
            os.system("mv mvAll.log.py .tmp.py")
        elif opt[ 0] =='--overwriteLog': logmode = 'overwrite'
        elif opt[ 0] =='--uniqueLog': logmode = 'unique'
        elif opt[ 0] =='--noLog': logmode = 'no'
        elif opt[ 0] =='--die': die = 1
        elif opt[ 0] =='--noSplash': noSplash = True
        elif opt[ 0] == '--customizer':
            customFile = opt[1]
        elif opt[ 0] == '--lib':
            libraries.append(opt[1])
        elif opt[ 0] in ('-i', '--interactive'):
            interactive = 1
        elif opt[ 0] in ('-d', '--dmode'):
            assert min([mo in dmodes for mo in opt[1].split('|')])==True
            dmode = opt[1]
        elif opt[ 0] in ('-c', '--cmode'):
            assert min([mo in cmodes for mo in opt[1].split('|')])==True
            cmode = opt[1]
        elif opt[ 0] in ('-v', '--vision'):
            if opt[1] in ('o', 'once'):
                visionarg = 'once'
            elif opt[1] in ('r', 'run'):
                visionarg = 'run'
        else:
            print "unknown option %s %s"%tuple(opt)
            print help_msg
            sys.exit( 1)
    
    text = 'Python executable     : '+sys.executable+'\n'
    if kw.has_key('AdtScriptPath'):
        text += 'ADT script                : '+ kw['AdtScriptPath'] +'\n'
    text += 'MGLTool packages '+'\n'
    
    
    from Support.path import path_text, release_path
    from Support.version import __version__
    from mglutil import __revision__
    version = __version__
    
    text += path_text
    text += version+': '+release_path
    
    path_data = text
    # if MGLPYTHONPATH environment variable exists - insert the specified path
    # into sys.path
    
    #if os.environ.has_key("MGLPYTHONPATH"):
    #    if sys.platform == "win32":
    #        mglPath = split(os.environ["MGLPYTHONPATH"], ";")
    #    else:
    #        mglPath = split(os.environ["MGLPYTHONPATH"], ":")
    #    mglPath.reverse()
    #    for p in mglPath:
    #        sys.path.insert(0, os.path.abspath(p))
            
        
    try:
        ##################################################################
        # Splash Screen
        ##################################################################
        print 'Run AutoDockTools from', __path__[0]    
        import Pmv
        image_dir = os.path.join( Pmv.__path__[0],'Icons','Images')
        copyright = """(c) 1999-2009 Molecular Graphics Laboratory, The Scripps Research Institute
    ALL RIGHTS RESERVED """
        authors = """Authors: Michel F. Sanner, Ruth Huey, Sargis Dallakyan,
    Sowjanya Karnati, William (Lindy) Lindstrom, Garrett M. Morris,
    Brian Norledge, Anna Omelchenko, Daniel Stoffler, Guillaume Vareille"""
        icon = os.path.join(Pmv.__path__[0],'Icons','64x64','adt.png')
        third_party = """Fast Isocontouring, Volume Rendering -- Chandrait Bajaj, UT Austin
    Adaptive Poisson Bolzman Solver (APBS) -- Nathan Baker Wash. Univ. St Louis
    GL extrusion Library (GLE) -- Linas Vepstas
    Secondary Structure Assignment (Stride) -- Patrick Argos EMBL
    Mesh Decimation (QSlim 2.0) -- Micheal Garland,  Univeristy of Illinois
    Tiled Rendering (TR 1.3) -- Brian Paul"""
        title = "AutoDockTools"
        #create a root and hide it
        from Tkinter import Tk
        root = Tk()
        root.withdraw()
    
        from mglutil.splashregister.splashscreen import SplashScreen
        from mglutil.splashregister.about import About
        about = About(image_dir=image_dir, third_party=third_party,
                      title=title, version=version, revision=__revision__,
                      path_data=path_data,
                      copyright=copyright, authors=authors, icon=icon)
        splash =  SplashScreen(about, noSplash=noSplash)
        from Pmv.moleculeViewer import MoleculeViewer

        mv = MoleculeViewer( logMode=logmode, customizer=customizer, master=root,
                             title=title, withShell= not interactive, verbose=False)
        mv.browseCommands('autotors4Commands', commands = None, 
                                        package = 'AutoDockTools')
        mv.browseCommands('autoflex4Commands', commands = None, 
                                     package = 'AutoDockTools')
        mv.browseCommands('autogpf4Commands', commands = None, 
                                        package = 'AutoDockTools')
        mv.browseCommands('autodpf4Commands', commands = None, 
                                        package = 'AutoDockTools')
        mv.browseCommands('autostart4Commands', commands = None, 
                                        package = 'AutoDockTools')
        mv.browseCommands('autoanalyze4Commands', commands = None, 
                                        package = 'AutoDockTools')
        mv.GUI.currentADTBar = 'AutoTools4Bar'
        setADTmode("AD4.0", mv)
        setADTmode("AD4.2", mv)
        mv.browseCommands('selectionCommands', package='Pmv')
        mv.GUI.naturalSize()
        mv.customize('_adtrc')
        
        mv.help_about = about
        
        font = mv.GUI.ROOT.option_get('font', '*')

        try:
            import Vision
            mv.browseCommands('visionCommands', commands=('vision',), topCommand=0)
            mv.browseCommands('coarseMolSurfaceCommands', topCommand=0)
            if hasattr(mv,'vision') and mv.vision.ed is None:
                mv.vision(log=0)
            else:
                # we address the global variable in vision
                Vision.ed = mv.vision.ed
        except ImportError:
            pass

        mv.GUI.ROOT.option_add('*font', font)  
        self = mv   
             
        #show the application after it built
        splash.finish()
        root.deiconify()
        globals().update(locals())

        mv.GUI.VIEWER.suspendRedraw = True
        cwd = os.getcwd() 
        #mv._cwd differs from cwd when 'Startup Directory' userpref is set
        os.chdir(mv._cwd)
        
        if dmode is not None or cmode is not None:
            # save current list of commands run when a molecule is loaded
            addCmds = mv.getOnAddObjectCmd()
        # remove them
            for c in addCmds:
                mv.removeOnAddObjectCmd(c[0])
            # set the mode
                setdmode(dmode, mv)
    
        if cmode is not None:
        # set the mode
            setcmode(cmode, mv)
    
        for a in args:
            if a[0]=='-':# skip all command line options
                continue
    
            elif (a[-10:]=='_pmvnet.py') or (a[-7:]=='_net.py'):  # Vision networks
                mv.browseCommands('visionCommands', commands=('vision',) )
                if mv.vision.ed is None:
                    mv.vision()
                mv.vision.ed.loadNetwork(a)
                if visionarg == 'run' or visionarg == 'once':
                    mv.vision.ed.currentNetwork.run()

            elif a[-3:]=='.py':     # command script
                print 'sourcing', a
                mv.source(a)
    
            elif a[-4:] in ['.pdb', '.pqr', 'pdbq', 'mol2', '.cif'] or a[-5:]=='pdbqs' or a[-5:]=='pdbqt':
                mv.readMolecule(a)
    
            else:
                print 'WARNING: unable to process %s command line argument'%a
        if again:
            mv.source(".tmp.py")
        if dmode is not None or cmode is not None:
            # get current list of commands run when a molecule is loaded
            cmds = mv.getOnAddObjectCmd()
            # remove them
            for c in cmds:
                mv.removeOnAddObjectCmd(c[0])
            # restore original list of commands
            for c in addCmds:
                apply( mv.addOnAddObjectCmd, c )
        
        mv.GUI.VIEWER.suspendRedraw = False
        os.chdir(cwd)
        if visionarg != 'once':
            if ownInterpreter is True:
                mod = __import__('__main__')
                mod.__dict__.update({'self':mv})
                if interactive:
                    sys.stdin = sys.__stdin__
                    sys.stdout = sys.__stdout__
                    sys.stderr = sys.__stderr__
                    import code
                    try: # hack to really exit code.interact 
                        code.interact( 'AutoDockTools Interactive Shell', local=mod.__dict__)
                    except:
                        pass
                elif not die:
                    #mv.GUI.pyshell.interp.locals = globals()
                    mv.GUI.pyshell.interp.locals = mod.__dict__
                    mv.GUI.ROOT.mainloop()
                mod.__dict__.pop('self')

    except:
        import traceback
        traceback.print_exc()
        raw_input("hit enter to continue")
        import sys
        sys.exit(1)


FlagCheck = 1
CRITICAL_DEPENDENCIES = ['MolKit', 'mglutil', 'Pmw', 'numpy', 'ViewerFramework', 'Pmv', 'DejaVu', 'opengltk']
NONCRITICAL_DEPENDENCIES = ['PyBabel','PyAutoDock','PIL', 'NetworkEditor', 'Vision','ZSI', 'Support', 'python_cluster']


