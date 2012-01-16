GLmol - Molecular Viewer on WebGL/Javascript

== About GLmol ==

GLmol is a 3D molecular viewer based on WebGL and Javascript. You can
embed molecular models in Web pages without using Java or
plugins. GLmol is open-source software licensed under LGPL3.

== Features ==

Currently GLmol has following features. More is coming...

* Read PDB file
* Read SDF/MOL file
* Load local files (File API is used; Safari is not supported)
* Load PDB files directly from RCSB PDB server 
* Rotate/Translate/Zoom model with mouse (touchpanel support is under development) 
* Fog & Slab
* Representations
    - Line (depiction of double/triple bonds is supported for SDF/MOL file)
    - Stick
    - Sphere(van der Waals radius or fixed radius)
    - Alpha carbon trace
    - Ribbon
    - Helix as Cylinder
    - Tube with radius reflecting B factor
    - Combination of above
* Special representations for Nucleic acid bases
    - Stick
    - Ladder
    - Line
* Coloring
    - By chain
    - By secondary structure (when defined in SHEET/HELIX records)
    - By Elements
    - Gradation (a.k.a chainbow)
    - Polar/Nonpolar
    - B factor
    - Custom
* Crystallography
    - Display unit cell
    - Show crystal packing (when defined in REMARK section)
* Display biological assembly (when defined in REMARK section)

== System requirements ==

GLmol runs on newer versions of Firefox, Chrome, Safari or
Opera. Internet Explorer is not supported because IE doesn't implement
WebGL. GLmol runs on Sony Ericsson's Android devices which support
WebGL. Support for Firefox Mobile is currently underway. Reportedly,
GLmol also runs on WebGL enabled safari in iOS.

== Troubleshooting ==

If you see only black screen and you are using

 Internet Explorer: sorry. IE doesn't support WebGL.
 Firefox (version 4 or later): try force enable WebGL. 
   https://wiki.mozilla.org/Blocklisting/Blocked_Graphics_Drivers#How_to_force-enable_blocked_graphics_features
 Chrome: try force enable WebGL.
   http://www.google.com/support/forum/p/Chrome/thread?tid=4b9244822aa2f2e0&hl=en
 Safari: enable WebGL.
   https://discussions.apple.com/thread/3300585?start=0&tstart=0

== How to embed ==

Currently, documentation is not ready.
Please examine "embedding-examplesEN.html"

Please note that Same-Origin-Policy applies to XmlHttpRequest so that
GLmol can load PDB files only from the same server as the program or
RCSB PDB. You can also embed whole PDB file in a HTML file.
 (see triiodotyrosine example)

== Contact ==

Project website is located at http://webglmol.sourceforge.jp/

Comments and suggestions are welcome at http://sourceforge.jp/projects/webglmol/forums/ or 
biochem_fan@users.sourceforge.jp 
