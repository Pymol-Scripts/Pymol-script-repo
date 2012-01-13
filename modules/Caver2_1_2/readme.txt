====== version  ========
Caver 2.1.2


====== requirements =======
Java Runtime Environment at least 1.6 (6) required 
http://java.sun.com/javase/downloads/


====== installation instructions: WINDOWS =======

1. Make sure java runtime environment (version at least 1.6) is installed on your system.
   If not sure, run command prompt and type "java -version" 

2. unzip "windows/setup.exe". Run and follow the instructions. Make sure you specify the directory with pyMol correctly when asked. Also make sure you specify your desired output directory correctly.

3. OPTIONAL (if you need to change default location of pyMol or Output dir)
   To specify default locations you can perform the following steps:
   
   3.0  Go to file: "PATH_TO_PYMOL/modules/pmg_tk/startup/Caver2_1_<version>.py" and edit it with text editor
   3.1  Replace  OUTPUT_LOCATION = "C:\\Caver2_1\\Output" (line 39) with your preferred output directory
   3.2  Replace  PYMOL_LOCATION = "C:\\Program Files\\DeLano Scientific\\PyMOL" (line 44) with your path to pyMol



====== installation instructions: LINUX & MAC =======

1. Make sure java runtime environment (version at least 1.6) is installed on your system.
   If not sure, run terminal and type "java -version" 

2. Extract "linux_mac/Caver2_1_<version>" folder to some destination.
   Example: /home/user/programs/Caver2_1_1/ 

3. Extract and edit "linux_mac/Caver2_1_<version>.py" with a text editor.
  3.1 replace "directory/where/jar/with/plugin/is/located" (line 46) with the path used in step 2.
      Note: check that Caver2_0.jar is located directly in the specified directory
  3.2 replace "/tmp" (line 41) with your preferred output directory

4. Run pyMol and click "Plugin -> Install plugin" and select the previously edited and saved "Caver2_1_<version>.py"

  
====== Changelog =======

Caver 2.1.2: tunnel spheres are ordered from the active site towards the surface
             minor bugfixes to unify resulting files (summary.txt, *.pqr, *.py)
Caver 2.1.1: added advanced input model support; improved installer