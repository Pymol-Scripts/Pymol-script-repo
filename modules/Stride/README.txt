Installing stride from source
http://webclu.bio.wzw.tum.de/stride/install.html

In order to install stride on a Unix/Linux system with gcc follow these steps. Other compilers need slight modifications of the Makefile.
Download the stride stand-alone program source code
Create a directory wherever/stride and copy the tar ball there
Unpack the tar ball: tar -zxf stride.tar.gz
Compile stride: make
Now you have a binary file stride in your directory, you can copy it to /usr/local/bin or add the current directory to your path to use it system-wide
