# meshscript
Visualize 3d meshes via scripting in scheme.

Content
-------
* [Introduction](#introduction)
* [Building](#building)
* [Examples](#examples)
     - [One bit render of 3D mesh to png](#one-bit-render-of-3d-mesh-to-png)
* [Glossary](#glossary)
* [Credits](#credits)

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
    
will load the mesh. The return value is an integer or id that represents this mesh. If you want to change any properties of the mesh, you'll need this id. Therefore it's probably better to load a mesh as

    (define id (load-mesh "D:/my_3d_files/rabbit.ply"))
    
so that you can query properties of this mesh. For instance

    (info id)
    
shows the number of vertices and triangles and other properties of the mesh.
If you want to visualize the mesh, you'll need to call

    (view-show!)
    
to open up the 3D view. The mesh will be displayed. You can hide the view again by closing it with the mouse or by calling

    (view-hide!)
    
in the REPL. If you want to keep the view, but hide your 3D model, you can call

    (hide! id)
    
In this case we've been writing our code directly in the REPL. You can also write your script in a separate file, and then call meshscript with this script file as argument. Meshscript will start, compile the script, run the script, and then start the REPL. All functions or defines from the script will be available in the REPL.    

Examples
--------

### One bit render of 3D mesh to png
![](images/porsche.png)
For this example I've used a free model of a Porsche 911 available via the link https://www.thingiverse.com/thing:2946486

    ; porsche 911 race car: https://www.thingiverse.com/thing:2946486
    
    ; load the mesh
    (define id (load-mesh "D:/stl/porsche/complete.stl"))
    
    ; set view properties
    (view-onebit-set! #t)
    (view-size-set! 1920 1024)
    (view-bg-set! 0 0 0)
    (view-shadow-set! #t)
    
    ; set coordinate system for the view
    (define cs '((0.593359 -0.220264 0.774216 11.876) (0.803842 0.212337 -0.555654 -2.86377) (-0.0420054 0.952049 0.303049 5.89276) (0 0 0 1)))
    (view-cs-set! cs)
    
    ; show the view window
    (view-show!)
    
    ; save the view as png image
    (view-export "D:/stl/porsche/porsche.png")





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
    	the y-axis, and by `z` degrees over the z-axis.
    
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
    	system. The coordinate system `cs` can be given as
    	a vector of size 16 in column major format or as a
    	list of lists in row major format.
    
    NAME
    	distance-map
    DESCRIPTION
    	(distance-map id1 id2 bool-signed) returns a list with
    	values that represent the distance between objects
    	with tag `id1` and `id2`. For each vertex of object
    	`id1` there is exactly one distance in the list. The
    	distance can be signed or unsigned, depending on the
    	boolean value that is given to `bool-signed`.
    
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
    	you want to use view-position in your script, as view-position
    	uses the data of the last render of the view.
    
    NAME
    	hide!
    DESCRIPTION
    	(hide! id) makes the object with tag `id` invisible.
    
    NAME
    	icp
    DESCRIPTION
    	(icp id1 id2 inlier-distance) returns the result of
    	the iterative closest point algorithm between objects
    	with tag `id1` and `id2`. This result is always a 4x4
    	transformation matrix. The iterative closest point
    	algorithm will only use correspondences between `id1`
    	and `id2` if their distance is smaller than `inlier-distance`.
    
    NAME
    	info
    DESCRIPTION
    	(info id) prints info on the object with tag `id`.
    
    NAME
    	jet
    DESCRIPTION
    	(jet lst) takes a list `lst` of values between 0 and
    	1 and returns a list of lists with (r g b) values.
    
    NAME
    	load-mesh
    DESCRIPTION
    	(load-mesh "stlfile.stl") loads the STL file and returns
    	an id. Similarly (load-mesh "objfile.obj") loads an
    	OBJ file and returns the id. Other input mesh formats
    	that are implemented are PLY and OFF.
    
    NAME
    	load-morphable-model
    DESCRIPTION
    	(load-morphable-model "model2019_fullHead.h5") loads
    	morphable models following the hdf5 file format as
    	used by the Basel Face Model project (https://faces.dmi.unibas.ch/bfm/bfm2019.html).
    	The other file format that can be read is meshscripts
    	own binary morphable model file format with extension
    	SSM.
    
    NAME
    	load-pointcloud
    DESCRIPTION
    	(load-pointcloud "pointcloud.ply") loads the PLY file
    	as point cloud and returns an id. Other file formats
    	are not yet supported.
    
    NAME
    	load-shape-predictor
    DESCRIPTION
    	(load-shape-predictor "filename") initializes the shape
    	predictor with the data given by "filename" and returns
    	the id. This is the dlib shape predictor (http://dlib.net).
    	The 68 points facial landmarks predictor data can be
    	downloaded from https://github.com/davisking/dlib-models
    
    NAME
    	make-mesh
    DESCRIPTION
    	(make-mesh vertices triangles) creates the mesh with
    	given `vertices` and `triangles`, and returns the id
    	of the created object. `vertices` should be a list
    	of lists of the form ((x y z) (x y z) ...) with x,y,z
    	floating point values, and `triangles` should be a
    	list of lists of the form ((a b c) (d e f) ...) with
    	a,b... fixnums referring to the vertex indices.
    
    NAME
    	make-pointcloud
    DESCRIPTION
    	(make-pointcloud vertices) creates the pointcloud with
    	given `vertices`, and returns the id of the created
    	object. `vertices` should be a list of lists of the
    	form ((x y z) (x y z) ...) with x,y,z floating point
    	values.
    
    NAME
    	marching-cubes
    DESCRIPTION
    	(marching-cubes bb dim isovalue fun) with `bb` representing
    	the bounding box of the form ((min_x max_x) (min_y
    	max_y) (min_z max_z)), `dim` representing the dimensions
    	of the form (width height depth), `isovalue` a flonum
    	representing the signed distance requested, and `fun`
    	representing the distance functions as a lambda function
    	accepting (x y z) values and returning a distance.
    
    NAME
    	matcap-set!
    DESCRIPTION
    	(matcap-set! id matcap-id) changes the matcap of the
    	object with tag `id`. The matcap is given by its id
    	matcap-id. Currently the matcaps in meshscript are
    	hardcoded. There are 4 available matcaps with ids 0,
    	1, 2, 3.
    
    NAME
    	mesh->pointcloud
    DESCRIPTION
    	(mesh->pointcloud id) converts the mesh with tag `id`
    	to a pointcloud.
    
    NAME
    	mesh-texture->vertexcolors
    DESCRIPTION
    	(mesh-texture->vertexcolors id) will return a list of
    	lists of the form ((r g b) (r g b) ... ). Each vertex
    	of the object with tag `id` has a corresponding (r
    	g b) value. This (r g b) value is obtained from the
    	texture of `id`, if available.
    
    NAME
    	morphable-model-coefficients-size
    DESCRIPTION
    	(morphable-model-coefficients-size mm_id) returns the
    	number of coefficients for the morphable model with
    	tag `mm_id`.
    
    NAME
    	morphable-model-shape-size
    DESCRIPTION
    	(morphable-model-shape_size mm_id) returns the shape
    	size for the morphable model with tag `mm_id`. This
    	is equal to the number of rows in the U matrix, where
    	a shape S is represented as S = mu + U*c, with mu the
    	average shape, and c the coefficients vector.
    
    NAME
    	morphable-model-sigma
    DESCRIPTION
    	(morphable-model-sigma mm_id idx) returns sigma for
    	the morphable model with tag `mm_id` at coefficient
    	index `idx`.
    
    NAME
    	morphable-model-coefficients
    DESCRIPTION
    	(morphable-model-coefficients mm_id) returns the list
    	of coefficients for the morphable model with tag `mm_id`.
    
    NAME
    	morphable-model-basic-shape-coefficients
    DESCRIPTION
    	(morphable-model-basic-shape-coefficients mm_id idx)
    	returns the list of coefficients of the `idx`-th shape
    	that was used to generate this morphable model. Not
    	all morphable models have this data. For instance the
    	Basel shape model does not contain this data.
    
    NAME
    	morphable-model-coefficients-set!
    DESCRIPTION
    	(morphable-model-coefficients-set! mm_id coeff) sets
    	the list of coefficients for the morphable model with
    	tag `mm_id`. Here `coeff` is a list of coefficient
    	values, and its size should equal (morphable-model-coefficients-size
    	mm_id).
    
    NAME
    	morphable-model->mesh
    DESCRIPTION
    	(morphable-model->mesh mm_id) converts the morphable
    	model with tag `mm_id` to a mesh and returns the new
    	id.
    
    NAME
    	morphable-model-color-coefficients-size
    DESCRIPTION
    	(morphable-model-color-coefficients-size mm_id) returns
    	the number of color coefficients for the morphable
    	model with tag `mm_id`.
    
    NAME
    	morphable-model-color-shape-size
    DESCRIPTION
    	(morphable-model-color-shape_size mm_id) returns the
    	shape size for the color part of the morphable model
    	with tag `mm_id`. This is equal to the number of rows
    	in the U_color matrix, where a the shape colors S_color
    	are represented as S_color = mu_color + U_color*c,
    	with mu_color the average shape colors, and c the color
    	coefficients vector.
    
    NAME
    	morphable-model-color-sigma
    DESCRIPTION
    	(morphable-model-color-sigma mm_id idx) returns sigma
    	for the colors of the morphable model with tag `mm_id`
    	at color coefficient index `idx`.
    
    NAME
    	morphable-model-color-coefficients
    DESCRIPTION
    	(morphable-model-color-coefficients mm_id) returns the
    	list of color coefficients for the morphable model
    	with tag `mm_id`.
    
    NAME
    	morphable-model-color-basic-shape-coefficients
    DESCRIPTION
    	(morphable-model-color-basic-shape-coefficients mm_id
    	idx) returns the list of color coefficients of the
    	`idx`-th color shape that was used to generate this
    	morphable model. Not all morphable models have this
    	data. For instance the Basel shape model does not contain
    	this data.
    
    NAME
    	morphable-model-color-coefficients-set!
    DESCRIPTION
    	(morphable-model-color-coefficients-set! mm_id coeff)
    	sets the list of color coefficients for the morphable
    	model with tag `mm_id`. Here `coeff` is a list of color
    	coefficient values, and its size should equal (morphable-model-color-coefficients-size
    	mm_id).
    
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
    	(npoint from to) computes the npoint-registration of
    	the set of 3d points in `from` to the set of 3d points
    	in `to`. The result is a 4x4 transformation matrix.
    	Here `from` and `to` are lists of lists of the form
    	((x y z) (x y z) ...) and `from` and `to` should have
    	the same amount of 3d points.
    
    NAME
    	pointcloud-normals-estimate!
    DESCRIPTION
    	(pointcloud-normals-estimate! id k) creates vertex normals
    	for the pointcloud with tag `id` by taking the `k`
    	nearest neighbours of each vertex, fit a least squares
    	plane through these points, and use the normal of this
    	plane as estimate.
    
    NAME
    	poisson
    DESCRIPTION
    	(poisson pc_id depth) applies Poisson surface reconstruction
    	to the pointcloud with tag `pc_id`. This is the screened
    	Poisson surface reconstruction algorithm by M. Kazhdan
    	and H. Hoppe. You have to provide the parameter `depth`
    	which represents the depth of the octree during Poisson
    	surface reconstruction.
    
    NAME
    	save
    DESCRIPTION
    	(save id "file.ext") writes the object with tag `id`
    	to file. The filetype is determined by the extension
    	that is given. You can export meshes to STL or PLY,
    	pointclouds to PLY, morphable models to SSM.
    
    NAME
    	shape-predict
    DESCRIPTION
    	(shape-predict sp_id (x y w h)) or (shape-predict sp_id
    	((x y w h) ...)) runs the shape predictor with tag
    	`sp_id` on the region defined by (x y w h) or on the
    	regions defined by ((x y w h) ...) in the current view
    	and returns the coordinates of the landmarks as a list
    	of lists. The predictor should be initialized with
    	load-shape-predictor.
    
    NAME
    	shape-predictor-horizontal-flip-set!
    DESCRIPTION
    	(shape-predictor-horizontal-flip-set! id #t/#f) toggles
    	horizontal flipping of the shape predictor given by
    	tag `id`.
    
    NAME
    	shape-predictor-link-to-ear-left-detector
    DESCRIPTION
    	(shape-predictor-link-to-ear-left-detector id) links
    	the shape predictor given by tag `id` to the ear left
    	detector. The result of this operation is that the
    	shape predictor is rendered automatically when the
    	left ear detector's automatic rendering is on. You
    	can turn on/off automatic rendering of the left ear
    	detector with the command (view-ear-left-detector-set!
    	#t/#f).
    
    NAME
    	shape-predictor-link-to-ear-right-detector
    DESCRIPTION
    	(shape-predictor-link-to-ear-right-detector id) links
    	the shape predictor given by tag `id` to the ear right
    	detector. The result of this operation is that the
    	shape predictor is rendered automatically when the
    	right ear detector's automatic rendering is on. You
    	can turn on/off automatic rendering of the right ear
    	detector with the command (view-ear-right-detector-set!
    	#t/#f).
    
    NAME
    	shape-predictor-link-to-face-detector
    DESCRIPTION
    	(shape-predictor-link-to-face-detector id) links the
    	shape predictor given by tag `id` to the face detector.
    	The result of this operation is that the shape predictor
    	is rendered automatically when the face detector's
    	automatic rendering is on. You can turn on/off automatic
    	rendering of the face detector with the command (view-face-detector-set!
    	#t/#f).
    
    NAME
    	shape-predictor-unlink
    DESCRIPTION
    	(shape-predictor-unlink id) unlinks the shape predictor
    	given by tag `id`, see shape-predictor-link-to-face-detector,
    	shape-predictor-link-to-ear-right-detector, or shape-predictor-link-to-ear-left-detector.
    
    NAME
    	show!
    DESCRIPTION
    	(show! id) makes the object with tag `id` visible.
    
    NAME
    	triangles
    DESCRIPTION
    	(triangles id) returns the triangles of object with
    	tag `id` as a list of lists of the form ((v0 v1 v2)
    	(v3 v4 v4) ...) where each sublist (v0 v1 v2) contain
    	the indices of the vertices that form a triangle. The
    	actual vertex positions can be obtained with the command
    	(vertices id).
    
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
    	(vertices id) returns the vertices of object with tag
    	`id` as a list of lists of the form ((x y z) (x y z)
    	...) where each sublist (x y z) is a 3d point representing
    	the position of that vertex.
    
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
    	view camera. The coordinate system `cs` can be given
    	as a vector of size 16 in column major format or as
    	a list of lists in row major format.
    
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
    	(view-position x y) returns the 3d position of the vertex
    	in the last render of the 3d view at coordinate (x,y).
    
    NAME
    	view-index
    DESCRIPTION
    	(view-index x y) returns the vertex index of the vertex
    	in the last render of the 3d view at coordinate (x,y).
    
    NAME
    	view-id
    DESCRIPTION
    	(view-id x y) returns the id of the object in the last
    	render of the 3d view at coordinate (x,y).
    
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
    	(exit) can be used in the input script to end meshscript,
    	so the REPL is skipped.
    
    
Credits
--------
- Image file formats: https://github.com/nothings/stb
- Iterative closest point: http://www.cvlibs.net/software/libicp/
- Marching cubes: http://paulbourke.net/geometry/polygonise/
- Morphable models: https://faces.dmi.unibas.ch/bfm/bfm2019.html
- Ply file import/export: http://w3.impa.br/~diego/software/rply/
- Poisson surface reconstruction: https://github.com/mkazhdan/PoissonRecon
- Shape detection and prediction: http://dlib.net/
- Software rendering (qbvh): https://www.embree.org/
