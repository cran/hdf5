hdf5 R package - Windows Read Me:
---------------------------------

To compile and then use the hdf5 library, you need
the Windows binary of the hdf5 library available at:
http://hdf.ncsa.uiuc.edu. We tested it with the
Windows binary version 1.4.3.

Download the zip package and unzip it somewhere
(for instance, in c:\temp; this dir will be called
<hdf5root> hereunder).


1) Compilation of the library:
------------------------------
- Unzip the current source files of the library in
<R dir>\src\hdf5, if it is not already done.

- Copy hdf5dll.dll from
<hdf5root>\5-143-win\c\release\dll to
<R dir>\src\hdf5\src

- Copy <hdf5root>\5-143-win\c\release\include to
<R dir>\src\hdf5\src\include (all files in this dir)

- Now you should compile successfully as usually with
Rcmd INSTALL %R_HOME%\src\hdf5

- You can zip the dir <R dir>\library\hdf5 for making
installable Windows binary package


2) Installation of the hdf5 library:
------------------------------------
You cannot run function in this package if you do not
put the corresponding hdf5dll.dll in your path...
WARNING: using a different version of the dll than
the one that was used for compilation could crash R!
