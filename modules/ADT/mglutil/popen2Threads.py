import sys, os
from threading import Thread
import Pmw, Tkinter
import Queue
try:
    import subprocess
except ImportError:
    import process as subprocess
    
if sys.platform=='win32':
    mswin = True
else:
    mswin = False

## TODO:
##     -kill button under windows
    
class SysCmdInThread(Thread):
    """This class is use the run a child process in a separate thread, so that
reading its stdout and stderr streams will no block the parent process.

The strings printed by the child to stdout and stderr are dispalyed in a Pmw
ScrolledText widget.  A kill button allows killing the job (after confirmation)
The 'Ok' button dissmisses the GUI. The 'Ok' button only becomes enabled after
the child process has completed.

WARNING: if the child process runs a python interpreter, pass the -u
command line argument to prevent python from buffering the output.

IMPLEMENTATION NOTES:
This uses the pyexpect module for posix systems and the popen5.process.py
reference implementation of PEP 324 PEP 324,  http://www.python.org/peps/pep-0324.html
The windows implementation requires the win32all extensions to be installed
http://sourceforge.net/project/showfiles.php?group_id=78018
"""
    def __init__(self, cmdstring, input=None, hasGui=True, shell = False):
        Thread.__init__(self, name="Producer")
        self.gui = hasGui
        self.input = input
        self.debug = 0 # to print debug message
        self.cmdstring = cmdstring
        self.output = Queue.Queue(-1) # Thread safe queue infinite size
        self.shell = shell
        if self.gui:
            self.root = Tkinter.Toplevel()
            if type(cmdstring) == list:
                cmdLine = ' '.join(cmdstring)
            else:
                cmdLine = cmdstring
            self.root.title(cmdLine)
            self.root.config(cursor='watch')
            fixedFont = Pmw.logicalfont('Fixed')
            self.stdoutTk = Pmw.ScrolledText( self.root, labelpos = 'nw',
                                    label_text='stdout', usehullsize = 1,
                                    hull_width = 700, hull_height = 500,
                                    text_wrap='none', text_font = fixedFont,
                                     text_foreground = 'blue', 
                                     text_background='white' )
            self.stdoutTk.pack(expand=True, fill='both')
            self.stdoutTk.insert('end', "Running: " + cmdLine +"\n")
            self.cursor = self.stdoutTk.component('text').config('cursor')[-1]
            self.stdoutTk.component('text').config(cursor='watch')
            f = Tkinter.Frame(self.root)
            f.pack(fill='x')            
            #Tkinter.Label(f, text='  ').pack(side='right') #spacer
            if not mswin:
                self.kill = Tkinter.Button(f, text='   Kill   ', 
                                           command=self.kill_cb)
                self.kill.pack(side='right')
            
            self.ok = Tkinter.Button(f, text='OK', command=self.ok_cb,
                                     state='disabled')
            self.ok.pack(side='right', expand=True, fill='x')

            self.update_me()
            
        
    def run(self):
        done = False
        try:
            self.com = subprocess.Popen(self.cmdstring, stdin=subprocess.PIPE,\
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT, shell=self.shell)
        except Exception, inst:
            # FIXME: Segmentation fault is not handled
            #print >>sys.stderr, "Execution failed:", inst
            self.output.put(inst)
            self.output.put(None)
            return
        inp, out, err = (self.com.stdin, self.com.stdout, self.com.stderr)
        # send data to stdin
        if self.input is not None:
            map( lambda x, f=inp: f.write("%s\n"%x), self.input)
            inp.close()
                
        while not done:
            data = out.readline()
            if data=='':
                done = True
            else:
                self.output.put(data)
                # only there for debug, you should access the stdout and stderr
                # by parsing self.output
                if self.debug:
                    print data
        self.output.put(None)

    def update_me(self):
        try:
            while 1:
                line = self.output.get_nowait()
                if line is None:
                    self.ok.configure(state='normal')
                    if not mswin:
                        self.kill.configure(state='disabled')
                    self.stdoutTk.component('text').config(cursor=self.cursor)
                    self.root.config(cursor='')
                    return
                else:
                    if type(line) == type('s'):
                        line = line.replace('\r','')
                    self.stdoutTk.insert('end', line)
                    self.stdoutTk.component('text').yview('end')      
                self.root.update_idletasks()
        except Queue.Empty:
            pass
        self.root.after(10, self.update_me)

    def kill_cb(self, event=None):
        if not mswin:
            import signal, Tkinter
            from SimpleDialog import SimpleDialog
            text = "Do you really want to kill this process?"
            d = SimpleDialog(Tkinter._default_root, text=text,
                             buttons=["Yes", "No"],
                             default=0, title="Kill process")
            result = d.go()
            if result==0:
                try:
                    os.kill(self.com.pid,signal.SIGKILL)
                except:
                    pass
        else:
            pass
        
    def ok_cb(self, event=None):
        self.root.destroy()
        #self.stdoutTk.destroy()
        #self.stderrTk.destroy()
        

#pgm = os.path.join('.', 'Binaries', sys.platform, 'rbox')

#cmd1 = SysCmdInThread('%s c d D2'%pgm)
#cmd1.start()
#cmd2 = SysCmdInThread('%s -u except.py'%sys.executable)
#cmd2.start()
#cmd3 = SysCmdInThread('%s -u bigLoop.py'%sys.executable)
#cmd3.start()
#cmd4 = SysCmdInThread('%s -u longJob.py'%sys.executable)#, hasGui=False)
#cmd4.start()
#cmd5 = SysCmdInThread('C:\\cygwin\\usr\\local\\bin\\msms -if test.xyzr -de 3.0')
#cmd5 = SysCmdInThread('msms -if 1gav.xyzr -de 3.0')
#cmd5.start()


#data = cmd1.output
#print "****************"
#for d in data: print d
#print "****************"

# WARNING HERE WE NEED to make sure cmd1. finished before we can run
#pgm = os.path.join('.', 'Binaries', sys.platform, 'qhull')
#cmd6 = SysCmdInThread('%s Qc s f Fx'%pgm, input=cmd1.output)
#cmd6.start()


## # QSLIM DECIMATION EXAMPLE
## f = open('1crn.smf')
## datain = f.readlines()
## f.close()
## from popen5 import process
## p = process.Popen('propslim 5000', 
##                   stdin=process.PIPE, stdout=process.PIPE)
## inp, out, err = (p.stdin, p.stdout, p.stderr)

## # send input
## map( lambda x, f=inp: f.write("%s\n"%x), datain)
## inp.close()

## data = out.readlines()

## vertices = []
## faces = []
## normals = []
## colors = []
## for l in data:
##     w = l.split()
##     if w[0]=='v':
##         vertices.append( [float(w[1]),float(w[2]),float(w[3])] )
##     elif w[0]=='f':
##         faces.append( [int(w[1])-1,int(w[2])-1,int(w[3])-1] )
##     if w[0]=='n':
##         normals.append( [float(w[1]),float(w[2]),float(w[3])] )
##     if w[0]=='c':
##         colors.append( [float(w[1]),float(w[2]),float(w[3])] )
## from DejaVu import Viewer
## vi = Viewer()
## from DejaVu.IndexedPolygons import IndexedPolygons
## pol = IndexedPolygons('decimated', vertices=vertices, faces=faces,
##                       vnormals=normals, materials=colors)
## vi.AddObject(pol)

