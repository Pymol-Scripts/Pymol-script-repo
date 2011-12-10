from idlelib import PyShell
#from idlelib.PyShell import PyShell, PyShellFileList, use_subprocess
from idlelib.EditorWindow import fixwordbreaks
import Tkinter, sys, os
def getShell(thread, rootTk = None, subprocess = False, debug=False,
             enable_shell=False, enable_edit=True):
    """
    This function creates and returns a shell PyShell instance
    required arguments:
    thread        --

    optional arguments:
    rootTk        --
    subprocess    -- boolean flag when set to True a the pyshell
                     runs a new interpreter in a sub process
                     when set to False it the user has access to the main
                     interpreter. (default = False)
    enable_shell  -- boolean flag when set to True a python shell is
                     created (by default True)
    enable_edit   -- boolean flag when set to True a edit shell is created
                     aware of the python syntax (by default False)
    debug         -- boolean flag when set to True starts the debugger when
                     the pyshell is created. (by default = False)
    """
    cmd = None
    script = None
    startup = False
    
#    try:
#        sys.ps1
#    except AttributeError:
#        sys.ps1 = '>>> '

    if hasattr(sys, 'ps1') is False:
        sys.ps1 = '>>> '

    global mainThread
    
    PyShell.use_subprocess = subprocess
    mainThread = thread
    for i in range(len(sys.path)):
        sys.path[i] = os.path.abspath(sys.path[i])

    pathx = []
    for dir in pathx:
        dir = os.path.abspath(dir)
        if not dir in sys.path:
            sys.path.insert(0, dir)

    global flist, root
    if rootTk is None: root = Tkinter.Tk()
    else: root = rootTk
    fixwordbreaks(root)

    flist = PyShell.PyShellFileList(root)
    if enable_edit:
        flist.new()
        if enable_shell :
            flist.open_shell()
    elif enable_shell:
        flist.pyshell = PyShell.PyShell()
        #flist.pyshell.begin()
    shell = flist.pyshell
    if debug:
       shell.open_debugger()
    return shell


