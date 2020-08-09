# meshscript
Visualize 3d meshes via scripting in scheme

Introduction
------------

Meshscript is a scriptable interface for visualizing and editing 3d meshes and point clouds.
![](images/meshscript.gif)

Building
--------

Meshscript depends on [Intel's TBB library](https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html), and on [SDL2](https://www.libsdl.org/download-2.0.php). Both TBB and SDL2 are not delivered with the code and thus need to be installed by the user.

On Windows you can download TBB's binaries from its website, and install them, preferably, in 
folder C:\Program Files\TBB. Another folder is also possible, but then you'll need to
adapt the CMakeLists.txt file and make it point to the correct location.
On Ubuntu you can simply run 
  sudo apt install libtbb-dev 
to install TBB.

To install SDL2 on Windows, download its sources from its website, and build with cmake. Next install SDL2 to folder C:\Program Files\SDL2. Again, another folder is fine, but then you'll need to adapt the CMakeLists.txt file. On Ubuntu run

    sudo apt-get install libsdl2-dev

to install SDL2.

Next a solution file / makefile can be generated with CMake. Use Visual Studio or make to build the code.

The scripting funcionality uses the [skiwi compiler](https://github.com/janm31415/skiwi), wich is a just in time scheme compiler. Some functionality still needs to be built by skiwi during startup. The code that skiwi needs resides in the scm subfolder of the libskiwi folder. The compiler expects that the environment variable SKIWI_MODULE_PATH exists and points to this folder. The scm folder can be placed anywhere on your harddrive as long as SKIWI_MODULE_PATH points to it. Make sure that you use only slashes (/) and not backslashes (\\) when entering the path in the SKIWI_MODULE_PATH variable. Also make sure that you end the path with a slash, e.g. "C:/skiwi/scm/".

As soon as SKIWI_MODULE_PATH is correctly initialised you are ready to go.

Glossary
--------

`(hide-view)` hides the 3d view

`(show-view)` shows the 3d view

`(load-mesh "stlfile.stl")` loads the stl file and returns an id. Similarly (load-mesh \"objfile.obj\") loads an obj file and return the id.

`(load-pc "pointcloud.ply")` loads the ply file as point cloud and returns an id.

`(set-cs id cs)` sets a new coordinate system for mesh id. The coordinate system cs can be given as a vector of size 16 in column major format or as a list of lists in row major format.

`(get-cs id)` returns the coordinate system for mesh id.

`(make-mesh vertices triangles)` plots the mesh with given vertices and triangles, and returns the id of the plotted object. vertices should be a list of lists of the form ((x y z) (x y z) ...) with x,y,z floating point values, and triangles should be a list of list of the form ((a b c) (d e f) ...) with a,b... fixnums referring to the vertex indices.

`(vertices-to-csv id "file.csv")` exports the vertices of mesh id to a csv file

`(triangles-to-csv id "file.csv")` exports the triangles of mesh id to a csv file

`(show id)` shows mesh id

`(hide id)` hides mesh id

`(set-color id matcap-id)` changes the matcap of mesh id. The matcap is given by its id matcap-id.

`(jet lst)` takes a list of values between 0 and 1 and returns a list of lists with (r g b) values

`(set-vertex-colors id clrlst)` sets vertex colors for mesh id. The vertex colors are given as a list of lists with (r g b) values.

`(rotate id x y z)` rotates mesh id by x degrees over the x-axis, by y degrees over the y-axis, and by z degrees over the z_axis

`(translate id x y z)` translates mesh id by vector (x y z)

`(set-edges #t/#f)` turns on/off rendering of edges

`(set-shading #t/#f)` turns on/off lighting

`(set-textured #t/#f)` turns on/off rendering of texture

`(set-shadow #t/#f)` turns on/off rendering of shadow

`(set-wireframe #t/#f)` turns on/off rendering of wireframe

`(set-one-bit #t/#f)` turns on/off one-bit rendering

`(set-image-size w h)` resizes the plotted image to size (w, h)

`(unzoom)` sets the camera to its initial position

`(export-image "image-file.png")` exports the current view to a png image

`(set-bg-color r g b)` changes the background color to (r g b).

`(get-position x y)` returns the 3D position of coordinate (x,y).
