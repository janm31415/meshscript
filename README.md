# meshscript
Visualize 3d meshes via scripting in scheme.

Introduction
------------

Meshscript is a scriptable interface for visualizing and editing 3d meshes and point clouds.
![](images/meshscript.gif)

Building
--------

Meshscript has two dependencies that are not delivered with this source code:
  - [Intel's TBB library](https://software.intel.com/content/www/us/en/develop/tools/threading-building-blocks.html)
  - [SDL2](https://www.libsdl.org/download-2.0.php)
  
While TBB is optional (see later), SDL2 is required.

##### SDL2
To install SDL2 on Windows, download its sources from its website, and build with CMake. Next install SDL2 to folder C:\Program Files\SDL2. Another folder is fine, but then you'll need to adapt this source codes CMakeLists.txt file as it assumes the location C:\Program Files\SDL2 for SDL2. 

On Ubuntu run

    sudo apt-get install libsdl2-dev

to install SDL2.

On MacOS download the SDL2 framework from the [SDL2](https://www.libsdl.org/) website and install in /Library/Frameworks/

##### TBB

TBB is only needed when the CMake variable JTK_THREADING is set to tbb. If you don't want to use tbb, you can set JTK_THREADING to std or ppl (Windows only).
However, if you want to use TBB, you'll need to install it first.

On Windows you can download TBB's binaries from its website, and install them, preferably, in 
folder C:\Program Files\TBB. Another folder is also possible, but then you'll need to
adapt the CMakeLists.txt file and make it point to the correct location.

On Ubuntu you can simply run 

    sudo apt install libtbb-dev 

to install TBB.

On MacOS you can run

    brew install tbb
    
If this gives an error in the sense of `Cannot write to /usr/local/Cellar` then you can solve this probably by updating your write privileges in this folder with the command `sudo chmod a+w /usr/local/Cellar`, and then try `brew` again.

##### Meshscript
Use CMake to create a Visual Studio solution file on Windows, makefile on Ubuntu, or XCode project on MacOS.


Glossary
--------

`(cs-ref id)` returns the coordinate system for the object with tag `id`.

`(cs-rotate! id x y z)` rotates the object with tag `id` by `x` degrees over the x-axis, by `y` degrees over the y-axis, and by `z` degrees over the z_axis.

`(cs-set! id cs)` sets a new coordinate system for the object with tag `id`. The coordinate system `cs` can be given as a vector of size 16 in column major format or as a list of lists in row major format.

`(cs-translate! id x y z)` translates the object with tag `id` by vector (x y z).

`(hide! id)` makes the object with tag `id` invisible.

`(jet lst)` takes a list of values between 0 and 1 and returns a list of lists with (r g b) values.

`(load-mesh "stlfile.stl")` loads the stl file and returns an id. Similarly (load-mesh \"objfile.obj\") loads an obj file and returns the id.

`(load-pointcloud "pointcloud.ply")` loads the ply file as point cloud and returns an id.

`(make-mesh vertices triangles)` plots the mesh with given vertices and triangles, and returns the id of the plotted object. Vertices should be a list of lists of the form ((x y z) (x y z) ...) with x,y,z floating point values, and triangles should be a list of lists of the form ((a b c) (d e f) ...) with a,b... fixnums referring to the vertex indices.

`(matcap-set! id matcap-id)` changes the matcap of the object with tag  `id`. The matcap is given by its id matcap-id.

`(show! id)` makes the object with tag `id` visible.

`(triangles->csv id "file.csv")` exports the triangles of the object with tag `id` to a csv file.

`(vertexcolors-set! id clrlst)` sets vertex colors for the object with tag `id`. The vertex colors are given as a list of lists with (r g b) values.

`(vertices->csv id "file.csv")` exports the vertices of the object with tag `id` to a csv file.

`(view-bg-set! r g b)` changes the background color to (r g b).

`(view-cs)` returns the coordinate system of the view camera.

`(view-cs-set! cs)` sets the coordinate system of the view camera.

`(view-edges-set! #t/#f)` turns on/off rendering of edges.

`(view-export "image-file.png")` exports the current view to a png image.

`(view-hide!)` hides the 3d view.

`(view-onebit-set! #t/#f)` turns on/off one-bit rendering.

`(view-ref x y)` returns the 3D position of coordinate (x,y).

`(view-shading-set! #t/#f)` turns on/off lighting.

`(view-shadow-set! #t/#f)` turns on/off rendering of shadow.

`(view-show!)` shows the 3d view.

`(view-size-set! w h)` resizes the plotted image to size (w, h).

`(view-textured-set! #t/#f)` turns on/off rendering of texture.

`(view-unzoom!)` sets the camera to its initial position.

`(view-wireframe-set! #t/#f)` turns on/off rendering of wireframe.
