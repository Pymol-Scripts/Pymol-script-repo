# This module handles License Agreements for add-on packages
# !!!CAUTION!!! improper use of this script can erase your hard drive
# $Header: /opt/cvs/python/packages/share1.5/mglutil/splashregister/license.py,v 1.12.2.1 2009/06/05 19:39:21 sargis Exp $
# $Id: license.py,v 1.12.2.1 2009/06/05 19:39:21 sargis Exp $
#

import Tkinter, Pmw, os, sys, tkMessageBox
import mglutil
from mglutil.util.misc import ensureFontCase
path = os.path.split(mglutil.__path__[0])[0]
def findFilePath(fileName, packageName):
    try:
        return findFilePathTry(fileName, packageName)
    except:
        return None
    
def findFilePathTry(fileName, packageName):
    """ Find the path to the file from the package packageName"""
    fullName = os.path.join(os.path.join(path,packageName), fileName)
    if os.path.exists(fullName):
        return fullName
    else:
        return None


tk_root = Tkinter.Tk()
tk_root.title("License Agreements for add-on packages")
label_1 = Tkinter.Label(text="License Agreements for add-on packages.",
                           bg ='pink',font=(ensureFontCase('Helvetica'), 16) )
label_1.pack(fill='both',expand=0)

sf = Pmw.ScrolledFrame(tk_root,hscrollmode='none',usehullsize=1,
                       hull_width = 500,hull_height = 470,
                       vertflex='expand',horizflex='elastic',)
sf.pack(fill='both',expand=1)
sframe = sf.interior()
balloon = Pmw.Balloon(tk_root)

acad_var = Tkinter.BooleanVar()

def toggle_acad_comm():
    global acad_var
    if acad_var.get():
        mslib_radio.invoke(0)
        stride_radio.invoke(0)
        UTPackages_radio.invoke(0)
    else:
        mslib_radio.invoke(1)
        stride_radio.invoke(1)
        UTPackages_radio.invoke(1)

n_f = Tkinter.Frame(sframe)
n_f.pack()
Tkinter.Label(n_f, font=(ensureFontCase('helvetica'), 14, 'bold'), text="""I you plan to use this software for:""" ).grid(row=0,column=0, sticky='w')

r1 = Tkinter.Radiobutton(n_f, text="Non commercial - Install all these packages", variable=acad_var, value=1, 
            command=toggle_acad_comm, state='active')
r1.grid(row=1,column=0, sticky='w')
Tkinter.Label(n_f, justify='left', fg = 'gray45', text = 
"""   Non commercial use of our software is unrestricted, however be sure 
   to cite the appropriate papers in publications. """
).grid(row=2,column=0, sticky='w')

r2 = Tkinter.Radiobutton(n_f, text="Commercial - Install none of these packages", variable=acad_var, value=0, 
            command=toggle_acad_comm, )
r2.grid(row=3,column=0, sticky='w')
r2.deselect()

Tkinter.Label(n_f, justify='left', fg = 'gray45',text = 
"""   The software components listed below have restrictions for commercial
   research. Please read and comply to the terms of the licenses of the  
   components.
""").grid(row=4,column=0, sticky='w')

mslib_agree = "Undecided"
found_mlsib = False
try:
    if findFilePath('','mslib'):
        found_mlsib = True
    else:
        found_mlsib = False  
        mslib_agree = True        
except ImportError:
    found_mlsib = False  
    mslib_agree = True
    
mslib_group = Pmw.Group(sframe, tag_pyclass = None)
if found_mlsib:
    mslib_group.pack(fill='both',expand=0)
label_mslib1 = Tkinter.Label(mslib_group.interior(), font=(ensureFontCase('helvetica'), 16),
                             text="MSLIB License Agreement", fg='red')

label_mslib1.pack(fill='both',expand=0)
balloon.bind(mslib_group, "MSLIB is used to compute molecular surfaces")
mslib_var =Tkinter.BooleanVar()
mslib_var.set(0)

def toggle_mslib():
    mslib_license_group.toggle()
    sf.reposition()
mslib_license_group = Pmw.Group(mslib_group.interior(),
                    tag_pyclass=Tkinter.Checkbutton,
                    tag_text='Show MSLIB License',
                    tag_command=toggle_mslib,
                    tag_variable=mslib_var)
if found_mlsib:
    mslib_license_group.pack(fill='both',expand=0)               
mslib_license_group.toggle()
mslib_text = """THE SCRIPPS RESEARCH INSTITUTE SOFTWARE LICENSE AGREEMENT
1. Grant Of Limited License: Software Use Restrictions. 
You are free to use the program (MSMS) 
ONLY FOR NON-COMMERCIAL SCIENTIFIC RESEARCH. 
This license is issued to you as an individual. 
For institution wide use or commercial use of the software you will 
have to follow the official registration procedure. 
To do so you can contact us by e-mail (sanner@scripps.edu), 
mail (Michel Sanner, TSRI,Molecular Biology Department, TCP26, 
10550 North Torrey Pines Road, 92037 La Jolla) or fax (858) 784-2341.

2. Credits: WE REQUEST THAT YOU AGREE TO ACKNOWLEDGE 
THE USE OF THE MSMS SOFTWARE THAT RESULTS IN 
ANY PUBLISHED WORK, INCLUDING SCIENTIFIC PAPERS, 
FILMS AND VIDEOTAPES BY CITING THE FOLLOWING REFERENCE:
Sanner, M.F., Spehner, J.-C., and Olson, A.J. (1996) 
Reduced surface: an efficient way to compute molecular surfaces. 
Biopolymers, Vol. 38, (3), 305-320.

3. Copying Restrictions: You will not sell or otherwise distribute commercially
these programs or derivatives to any other party, whether with or without
consideration.

4. Ownership of Software: You will not obtain, and will not attempt to 
obtain copyright coverage thereon without the express purpose written 
consent of The Scripps Research Institute.

5. Limitation Of Liability: You will hold harmless from 
all or any expenses you may incur as a result of or arising from your use, 
direct or indirect, of these materials. You understand that no other right 
or license to this program or its derivatives is granted or implied as 
a result of our transmission of same to you. 
"""

def mslib_callback(tag):
    global mslib_agree
    if tag == "Agree to license and install":
        mslib_agree = True
    else:
        mslib_agree = False
    
label_mslib2 = Tkinter.Label(mslib_license_group.interior(),text=mslib_text)
if found_mlsib:
    label_mslib2.pack(fill='both',expand=0,side='left')
mslib_radio = Pmw.RadioSelect(mslib_group.interior(),
                command = mslib_callback,
                labelpos = 'e',
                frame_borderwidth = 1,
                frame_relief = 'ridge',
        )
mslib_radio.add("Agree to license and install")
mslib_radio.add("Disagree and do not install")
if found_mlsib:
    mslib_radio.pack(fill='both',expand=0)

stride_agree = "Undecided"
found_stride = False

try:
    if findFilePath('','stride'):
        found_stride = True
    else:
        found_stride = False  
        stride_agree = True        
except ImportError:
    found_stride = False  
    stride_agree = True

stride_group = Pmw.Group(sframe, tag_pyclass = None)
if found_stride:
    stride_group.pack(fill='both',expand=0)
label_stride1 = Tkinter.Label(stride_group.interior(), 
text="STRIDE License Agreement", fg='red',font=("Helvetica", 16))
if found_stride:
    label_stride1.pack(fill='both',expand=0)
balloon.bind(stride_group, "STRIDE is used to compute secondary structures.")
stride_var =Tkinter.BooleanVar()
stride_var.set(0)
def toggle_stride():
    stride_license_group.toggle()
    sf.reposition()
stride_license_group = Pmw.Group(stride_group.interior(),
                    tag_pyclass=Tkinter.Checkbutton,
                    tag_text='Show STRIDE license',
                    tag_command=toggle_stride,
                    tag_variable=stride_var)
if found_stride:
    stride_license_group.pack(fill='both',expand=0)               
stride_license_group.toggle()
stride_text = """STRIDE: Protein secondary structure assignment from atomic coordinates
Dmitrij Frishman & Patrick Argos
All rights reserved, whether the whole or part of the program is concerned. 
Permission to use, copy, and modify this software and its documentation 
is granted for academic use, provided that:

1. This copyright notice appears in all copies of the software 
and related documentation.

2. The reference given below (Frishman and Argos, 1995) must be cited 
in any publication of scientific results based in part or completely 
on the use of the program.

3. Bugs will be reported to the authors. The use of the software in 
commercial activities is not allowed without a prior written commercial 
license agreement.

WARNING:
STRIDE is provided 'as-is' and without warranty of any kind, express, 
implied or otherwise, including without limitation any warranty of 
merchantability or fitness for a particular purpose. 
In no event will the authors be liable for any special, incidental, 
indirect or consequential damages of any kind, or any damages whatsoever 
resulting from loss of data or profits, whether or not advised of the 
possibility of damage, and on any theory of liability, arising out of 
or in connection with the use or performance of this software.

For calculation of the residue solvent accessible area the program 
NSC is used and was kindly provided by 
Dr. F.Eisenhaber (EISENHABER@EMBL-HEIDELBERG.DE). 
Please direct to him all questions concerning specifically 
accessibility calculations.

Reference:
Frishman,D & Argos,P. (1995) Knowledge-based secondary 
structure assignment. Proteins: structure, function and genetics, 
23, 566-579. """

label_stride2 = Tkinter.Label(stride_license_group.interior(),text=stride_text)
if found_stride:
    label_stride2.pack(fill='both',expand=0,side='left')

def stride_callback(tag):
    global stride_agree
    if tag == "Agree to license and install":
        stride_agree = True
    else:
        stride_agree = False

stride_radio = Pmw.RadioSelect(stride_group.interior(),
                labelpos = 'e',
                command = stride_callback,
                frame_borderwidth = 1,
                frame_relief = 'ridge',
        )
stride_radio.add("Agree to license and install")
stride_radio.add("Disagree and do not install")
if found_stride:
    stride_radio.pack(fill='both',expand=0)
#stride_radio.button(0).config(bg="pink")

ut_agree = "Undecided"
found_UTpackages = False

try:
    if findFilePath('','UTpackages'):
        found_UTpackages = True
    else:
        found_UTpackages = False  
        ut_agree = True        
except ImportError:
    found_UTpackages = False  
    ut_agree = True

UTPackages_group = Pmw.Group(sframe, tag_pyclass = None)
if found_UTpackages:
    UTPackages_group.pack(fill='both',expand=0)

label_UTPackages1 = Tkinter.Label(UTPackages_group.interior(), fg='red', 
                  font=("Helvetica", 16), text="UTPackages License Agreement")
if found_UTpackages:
    label_UTPackages1.pack(fill='both',expand=0)
balloon.bind(UTPackages_group, "UTPackages is used for isocontouring and volume rendering")
UTPackages_var =Tkinter.BooleanVar()
UTPackages_var.set(0)
def toggle_UTPackages():
    UTPackages_license_group.toggle()
    sf.reposition()
UTPackages_license_group = Pmw.Group(UTPackages_group.interior(),
                    tag_pyclass=Tkinter.Checkbutton,
                    tag_text='Show UTPackages license',
                    tag_command=toggle_UTPackages,
                    tag_variable=UTPackages_var)
if found_UTpackages:
    UTPackages_license_group.pack(fill='both',expand=0)               
UTPackages_license_group.toggle()
UTPackages_text = """TERMS AND CONDITIONS FOR COPYING, 
DISTRIBUTION AND MODIFICATION
1. This software is the copyright of 
THE UNIVERSITY OF TEXAS AT AUSTIN, 2005. 

2. The software is available under multiple licenses.

3. For non commercial educational and non commercial academic use, 
the software including source code, interface definitions and compile 
scripts are freely available. Any distribution of code, library 
or executables which contains any modules from this software should 
contain, conspicously and appropriately, on each copy, this copyright 
and license notice and should be freely distributed.

4. If you wish to incorporate parts of the Library into other free programs 
whose distribution conditions are incompatible with these, write to the 
author to ask for permission. For software which is copyrighted by the 
Free Software Foundation, write to the Free Software Foundation. 

5. For any other purpose, including commercial purposes, 
please contact The University of  Texas at Austin for a different license.

6. Credits: 

This software has been developed at the Computational and 
Visualization center under 
    
    Dr Chandrajit Bajaj,
    Computational Applied Mathematics Chair in Visualization,
    Professor of Computer Sciences,
    Director of Center for Computational Visualization,
    Department of Computer Sciences & 
    The Institute of Computational Engineering and Sciences,
     Center for Computational Visualization,
    201 East 24th Street, ACES 2.324A,
    1 University Station, C0200,
    Austin, TX 78712-0027.

We request that you agree to acknowledge the use of the software that 
results in any published work, including scientific papers, films and 
videotapes by citing the references listed below

For UTblur:

C. Bajaj, V. Siddavanahalli
Fast Feature Adaptive Surfaces and Derivatives Computation for 
Volumetric Particle Data
ICES and CS technical reports, The University of Texas at Austin, 2005.

For UTisocontour:

C. Bajaj, V. Pascucci, D. Schikore
Accelerated IsoContouring of Scalar Fields
Data Visualization Techniques, edited by C. Bajaj, John Wiley and Sons (1998).

C. Bajaj, V. Pascucci, D. Schikore 
Fast Isocontouring for Improved Interactivity 
Proceedings: ACM Siggraph/IEEE Symposium on Volume Visualization, 
ACM Press, (1996), San Francisco, CA. Pages: 39-46 99

For UTmesh:

Y. Zhang, C. Bajaj
Adaptive and Quality Quadrilateral/Hexahedral Meshing from 
Volumetric Data Computer Methods in Applied Mechanics 
and Engineering (CMAME), in press, 2005.

Y. Zhang, C. Bajaj, B-S. Sohn
3D Finite Element Meshing from Imaging Data
The special issue of Computer Methods in Applied Mechanics and 
Engineering (CMAME) on Unstructured
Mesh Generation, 194(48-49):5083-5106, 2005.

For UTsdf:

C. Bajaj, V. Siddavanahalli
An Adaptive Grid Based Method for Computing Molecular Surfaces 
and Properties ICES and CS Technical Reports, 
The University of Texas at Austin, 2005.

For UTvolrend:

C. Bajaj, Z. Yu, M. Auer 
Volumetric Feature Extraction and Visualization of Tomographic 
Molecular Imaging. Journal of Structural Biology, 
Volume 144, Issues 1-2, October 2003, Pages 132-143

For UTmolderivatives:

C. Bajaj, V. Siddavanahalli
Fast Feature Adaptive Surfaces and Derivatives Computation for Volumetric 
Particle Data ICES and CS Technical Reports, 
The University of Texas at Austin, 2005.

7. No warranty 

7a. Because the library is licensed free of charge, there is no warranty 
for the library, to the extent permitted by applicable law. except when 
otherwise stated in writing the copyright holders and/or other parties 
provide the library "as is" without warranty of any kind, either 
expressed or implied, including, but not limited to, the implied warranties 
of merchantability and fitness for a particular purpose. the entire risk as 
to the quality and performance of the library is with you. should the library 
prove defective, you assume the cost of all necessary servicing,
repair or correction. 

7b. In no event unless required by applicable law or agreed to in writing 
will any copyright holder, or any other party who may modify and/or 
redistribute the library as permitted above, be liable to you for damages, 
including any general, special, incidental or consequential damages 
arising out of the use or inability to use the library (including but 
not limited to loss of data or data being rendered inaccurate or 
losses sustained by you or third parties or a failure of the 
library to operate with any other software), even if such holder 
or other party has been advised of the possibility of such damages.  """

label_UTPackages2 = Tkinter.Label(UTPackages_license_group.interior(),text=UTPackages_text)
if found_UTpackages:
    label_UTPackages2.pack(fill='both',expand=0,side='left')


def ut_callback(tag):
    global ut_agree
    if tag == "Agree to license and install":
        ut_agree = True
    else:
        ut_agree = False

UTPackages_radio = Pmw.RadioSelect(UTPackages_group.interior(),
                labelpos = 'e',
                command = ut_callback,
                frame_borderwidth = 1,
                frame_relief = 'ridge')
UTPackages_radio.add("Agree to license and install")
UTPackages_radio.add("Disagree and do not install")
if found_UTpackages:
    UTPackages_radio.pack(fill='both',expand=0)
#UTPackages_radio.button(0).config(bg="pink")
def finish():
    global mslib_agree, stride_agree, ut_agree
    check_l = True
    txt = "Please click on Agree or Disagree button(s) next to: \n\n"
    if mslib_agree == "Undecided":
        txt += "    MSLIB License Agreement\n"
        check_l = False
    if stride_agree == "Undecided":
        txt += "    STRIDE License Agreement\n"
        check_l = False
    if ut_agree == "Undecided":
        txt += "    UTPackages License Agreement\n"
        check_l = False        
    if not check_l:
        tkMessageBox.showinfo("Choose components to install", txt)
        return

    if not mslib_agree:
            mslib = findFilePath('','mslib')
            if mslib:
                for root, dirs, files in os.walk(mslib, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                        for name in dirs:
                            try:
                                os.rmdir(os.path.join(root, name))
                            except:
                                pass
                try:
                    os.rmdir(mslib)
                except:
                    pass
            
    if not stride_agree:
            stride = findFilePath('','stride')
            if stride:
                for root, dirs, files in os.walk(stride, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                        for name in dirs:
                            try:
                                os.rmdir(os.path.join(root, name))
                            except:
                                pass
                try:
                    os.rmdir(stride)
                except:
                    pass

    if not ut_agree:
            ut = findFilePath('','UTpackages')
            if ut:
                for root, dirs, files in os.walk(ut, topdown=False):
                    for name in files:
                        os.remove(os.path.join(root, name))
                        for name in dirs:
                            try:
                                os.rmdir(os.path.join(root, name))
                            except:
                                pass
                try:
                    os.rmdir(ut)
                except:
                    pass
    tk_root.destroy()

s_label = Tkinter.Label(tk_root,
                        text='            Please make your selection and click')
s_label.pack(side='left')
Close_button = Tkinter.Button(tk_root,text='    Continue   ', command=finish)
Close_button.pack(side='left')

tk_root.protocol("WM_DELETE_WINDOW", finish)
tk_root.mainloop()
