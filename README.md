# meshscript
Visualize 3d meshes via scripting in scheme.

Introduction
------------

Meshscript is a scriptable interface for visualizing and editing 3d meshes and point clouds.
![](images/meshscript.gif)

Building
--------
First of all, meshscript uses submodules, so don't forget to also call

     git submodule update --init

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

TBB is only needed when the CMake variable JTK_THREADING is set to tbb. If you don't want to use tbb, you can set JTK_THREADING to std or ppl (ppl will not work on Ubuntu or MacOS).
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

Basics
------

When you run meshscript you'll find yourself in a scheme REPL. I'm using my own [skiwi](https://github.com/janm31415/skiwi) compiler, which is a JIT compiler that translates scheme code to machine code and then executes the code. The scheme compiler is R4RS-compliant and almost R5RS-compliant.

Apart from regular scheme syntax, you can also call meshscript related functions. They are described in more detail in the glossary below. By typing

    ,external
    
in the scheme REPL you'll get an overview of all the meshscript related functionality. If you are looking for information on a function that starts with `load`, you can type

    ,external load
    
in the REPL to find all meshscript functions that start with `load`.

Let's start with an example. Suppose you want to visualize a PLY-file, then the following call

    (load-mesh "D:/my_3d_files/rabbit.ply")
    
will load the mesh. The return value is an integer or id that represents this mesh from now. If you want to change any properties of this mesh, you'll need this id. Therefore it's probably better to load a mesh as

    (define id (load-mesh "D:/my_3d_files/rabbit.ply"))
    
so that you can query properties of this mesh, e.g.

    (info id)
    
If you want to visualize this mesh, you'll need to call

    (view-show!)
    
to open up the 3D view. The mesh will be displayed. You can hide the view again by closing it with the mouse or by calling

    (view-hide!)
    
in the REPL. If you want to keep the view, but hide your 3D model, you can call

    (hide! id)
    
In this case we've been writing our code directly in the REPL. You can also write your script in a separate file, and then call meshscript with this script file as argument. Meshscript will start, compile the script, run the script, and then start the REPL. All functions or defines will be available in the REPL to further investigate.    

Glossary
--------

Below follows a dump of all the meshscript methods so far. 

    NAME
    	cs
    DESCRIPTION
    	(cs id) returns the coordinate system for the object
    	with tag `id`.
    
    NAME
    	cs-apply!
    DESCRIPTION
    	(cs-apply! id) transforms the vertices of object with
    	tag `id` by its coordinate system, and sets its coordinate
    	system to the world.
    
    NAME
    	cs-rotate!
    DESCRIPTION
    	(cs-rotate! id x y z) rotates the object with tag `id`
    	by `x` degrees over the x-axis, by `y` degrees over
    	the y-axis, and by `z` degrees over the z_axis.
    
    NAME
    	cs-set!
    DESCRIPTION
    	(cs-set! id cs) sets a new coordinate system for the
    	object with tag `id`. The coordinate system `cs` can
    	be given as a vector of size 16 in column major format
    	or as a list of lists in row major format.
    
    NAME
    	cs-translate!
    DESCRIPTION
    	(cs-translate! id x y z) translates the object with
    	tag `id` by vector (x y z).
    
    NAME
    	cs-premultiply!
    DESCRIPTION
    	(cs-premultiply! id cs) premultiplies the coordinate
    	system of the object with tag `id` by the input coordinate
    	system.
    
    NAME
    	distance-map
    DESCRIPTION
    	(distance-map id1 id2 bool-signed)
    
    NAME
    	ear-right-detect
    DESCRIPTION
    	(ear-right-detect) runs the ear detector on the current
    	view and returns a list of lists of the form ((x y
    	w h) ...) where (x y w h) represents a rectangle containing
    	the right ear starting in corner (x,y) and with sizes
    	(w,h).
    
    NAME
    	ear-left-detect
    DESCRIPTION
    	(ear-left-detect) runs the ear detector on the current
    	view and returns a list of lists of the form ((x y
    	w h) ...) where (x y w h) represents a rectangle containing
    	the left ear starting in corner (x,y) and with sizes
    	(w,h).
    
    NAME
    	face-detect
    DESCRIPTION
    	(face-detect) runs the face detector on the current
    	view and returns a list of lists of the form ((x y
    	w h) ...) where (x y w h) represents a rectangle containing
    	the face starting in corner (x,y) and with sizes (w,h).
    
    NAME
    	force-redraw
    DESCRIPTION
    	(force-redraw) redraws the canvas. This is useful if
    	you want to use view-position.
    
    NAME
    	hide!
    DESCRIPTION
    	(hide! id) makes the object with tag `id` invisible.
    
    NAME
    	icp
    DESCRIPTION
    	(icp id1 id2 inlier-distance) returns the result of
    	iterative closest point as coordinate system.
    
    NAME
    	info
    DESCRIPTION
    	(info id) prints info on the object with tag `id`.
    
    NAME
    	jet
    DESCRIPTION
    	(jet lst) takes a list of values between 0 and 1 and
    	returns a list of lists with (r g b) values.
    
    NAME
    	load-mesh
    DESCRIPTION
    	(load-mesh "stlfile.stl") loads the stl file and returns
    	an id. Similarly (load-mesh "objfile.obj") loads an
    	obj file and returns the id.
    
    NAME
    	load-morphable-model
    DESCRIPTION
    	(load-mesh "stlfile.stl") loads the stl file and returns
    	an id. Similarly (load-mesh "objfile.obj") loads an
    	obj file and returns the id.
    
    NAME
    	load-pointcloud
    DESCRIPTION
    	(load-pointcloud "pointcloud.ply") loads the ply file
    	as point cloud and returns an id.
    
    NAME
    	load-shape-predictor
    DESCRIPTION
    	(load-shape-predictor "filename") initializes the shape
    	predictor with the data given by "filename" and returns
    	the id.
    
    NAME
    	make-mesh
    DESCRIPTION
    	(make-mesh vertices triangles) plots the mesh with given
    	vertices and triangles, and returns the id of the plotted
    	object. Vertices should be a list of lists of the form
    	((x y z) (x y z) ...) with x,y,z floating point values,
    	and triangles should be a list of lists of the form
    	((a b c) (d e f) ...) with a,b... fixnums referring
    	to the vertex indices.
    
    NAME
    	marching-cubes
    DESCRIPTION
    	(marching-cubes bb dim isovalue fun) with bb of the
    	form ((min_x max_x) (min_y max_y) (min_z max_z)), dim
    	of the form (width height depth), isovalue a flonum,
    	fun a lambda function accepting (x y z) values and
    	returning a distance.
    
    NAME
    	matcap-set!
    DESCRIPTION
    	(matcap-set! id matcap-id) changes the matcap of the
    	object with tag `id`. The matcap is given by its id
    	matcap-id.
    
    NAME
    	mesh-texture->vertexcolors
    NAME
    	morphable-model-coefficients-size
    NAME
    	morphable-model-shape-size
    NAME
    	morphable-model-sigma
    NAME
    	morphable-model-coefficients
    NAME
    	morphable-model-basic-shape-coefficients
    NAME
    	morphable-model-coefficients-set!
    NAME
    	morphable-model->mesh
    NAME
    	morphable-model-color-coefficients-size
    NAME
    	morphable-model-color-shape-size
    NAME
    	morphable-model-color-sigma
    NAME
    	morphable-model-color-coefficients
    NAME
    	morphable-model-color-basic-shape-coefficients
    NAME
    	morphable-model-color-coefficients-set!
    NAME
    	morphable-model-fit-indices!
    DESCRIPTION
    	(morphable-model-fit-indices! mm_id indices positions)
    
    NAME
    	morphable-model-fit!
    DESCRIPTION
    	(morphable-model-fit! mm_id mesh_id)
    
    NAME
    	npoint
    DESCRIPTION
    	npoint
    
    NAME
    	poisson
    DESCRIPTION
    	(poisson pc_id depth)
    
    NAME
    	save
    DESCRIPTION
    	(save id "file.ext")
    
    NAME
    	shape-predict
    DESCRIPTION
    	(shape-predict sp_id (x y w h)) or (shape-predict sp_id
    	((x y w h) ...)) runs the shape predictor with tag
    	sp_id on the region defined by (x y w h) or on the
    	regions defined by ((x y w h) ...) in the current view
    	and returns the coordinates of the landmarks as a list
    	of lists. The predictor should be initialized with
    	load-shape-predictor.
    
    NAME
    	shape-predictor-horizontal-flip-set!
    DESCRIPTION
    	(shape-predictor-horizontal-flip-set! id #t/#f) toggles
    	horizontal flipping of the shape predictor given by
    	tag id.
    
    NAME
    	shape-predictor-link-to-face-detector
    DESCRIPTION
    	(shape-predictor-link-to-face-detector id) links the
    	shape predictor given by tag id to the face detector
    
    NAME
    	shape-predictor-link-to-ear-right-detector
    DESCRIPTION
    	(shape-predictor-link-to-ear-right-detector id) links
    	the shape predictor given by tag id to the ear right
    	detector
    
    NAME
    	shape-predictor-link-to-ear-left-detector
    DESCRIPTION
    	(shape-predictor-link-to-ear-left-detector id) links
    	the shape predictor given by tag id to the ear left
    	detector
    
    NAME
    	shape-predictor-unlink
    DESCRIPTION
    	(shape-predictor-unlink id) unlinks the shape predictor
    	given by tag id.
    
    NAME
    	show!
    DESCRIPTION
    	(show! id) makes the object with tag `id` visible.
    
    NAME
    	triangles
    DESCRIPTION
    	(triangles id)
    
    NAME
    	triangles->csv
    DESCRIPTION
    	(triangles->csv id "file.csv") exports the triangles
    	of the object with tag `id` to a csv file.
    
    NAME
    	vertexcolors-set!
    DESCRIPTION
    	(vertexcolors-set! id clrlst) sets vertex colors for
    	the object with tag `id`. The vertex colors are given
    	as a list of lists with (r g b) values.
    
    NAME
    	vertices
    DESCRIPTION
    	(vertices id)
    
    NAME
    	vertices->csv
    DESCRIPTION
    	(vertices->csv id "file.csv") exports the vertices of
    	the object with tag `id` to a csv file.
    
    NAME
    	view-bg-set!
    DESCRIPTION
    	(view-bg-set! r g b) changes the background color to
    	(r g b).
    
    NAME
    	view-cs
    DESCRIPTION
    	(view-cs) returns the coordinate system of the view
    	camera.
    
    NAME
    	view-cs-set!
    DESCRIPTION
    	(view-cs-set! cs) sets the coordinate system of the
    	view camera.
    
    NAME
    	view-edges-set!
    DESCRIPTION
    	(view-edges-set! #t/#f) turns on/off rendering of edges.
    
    NAME
    	view-export
    DESCRIPTION
    	(view-export "image-file.png") exports the current view
    	to a png image.
    
    NAME
    	view-ear-left-detector-set!
    DESCRIPTION
    	(view-ear-left-detector-set! #t/#f) turns on/off rendering
    	of the left ear detector result.
    
    NAME
    	view-ear-right-detector-set!
    DESCRIPTION
    	(view-ear-right-detector-set! #t/#f) turns on/off rendering
    	of the right ear detector result.
    
    NAME
    	view-face-detector-set!
    DESCRIPTION
    	(view-face-detector-set! #t/#f) turns on/off rendering
    	of the face detector result.
    
    NAME
    	view-hide!
    DESCRIPTION
    	(view-hide!) hides the 3d view.
    
    NAME
    	view-onebit-set!
    DESCRIPTION
    	(view-onebit-set! #t/#f) turns on/off one-bit rendering.
    
    NAME
    	view-position
    DESCRIPTION
    	(view-position x y) returns the 3D position of coordinate
    	(x,y).
    
    NAME
    	view-index
    DESCRIPTION
    	(view-index x y) returns the vertex index of coordinate
    	(x,y).
    
    NAME
    	view-id
    DESCRIPTION
    	(view-id x y) returns the id of the object at coordinate
    	(x,y).
    
    NAME
    	view-shading-set!
    DESCRIPTION
    	(view-shading-set! #t/#f) turns on/off lighting.
    
    NAME
    	view-shadow-set!
    DESCRIPTION
    	(view-shadow-set! #t/#f) turns on/off rendering of shadow.
    
    NAME
    	view-show!
    DESCRIPTION
    	(view-show!) shows the 3d view.
    
    NAME
    	view-size-set!
    DESCRIPTION
    	(view-size-set! w h) resizes the plotted image to size
    	(w, h).
    
    NAME
    	view-textured-set!
    DESCRIPTION
    	(view-textured-set! #t/#f) turns on/off rendering of
    	texture.
    
    NAME
    	view-unzoom!
    DESCRIPTION
    	(view-unzoom!) sets the camera to its initial position.
    
    NAME
    	view-wireframe-set!
    DESCRIPTION
    	(view-wireframe-set! #t/#f) turns on/off rendering of
    	wireframe.
    
    NAME
    	exit
    DESCRIPTION
    	(exit) can be used in the input script to end meshscript.
    
 
