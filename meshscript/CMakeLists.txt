if (WIN32)
  set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include" CACHE PATH "the include path for TBB")
  set(TBB_LIBRARIES "C:/Program Files/TBB/lib/intel64/vc14/tbb.lib" "C:/Program Files/TBB/lib/intel64/vc14/tbb_debug.lib" CACHE FILEPATH "the tbb lib files")
  set(OPENGL_LIBRARIES opengl32.lib CACHE FILEPATH "opengl lib file")
endif (WIN32)

if (UNIX)
  if (APPLE)
    set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include" CACHE PATH "the include path for TBB")
    set(TBB_LIBRARIES "/usr/local/lib/libtbb.dylib" CACHE FILEPATH "the tbb .dylib file")
    set(OPENGL_LIBRARIES "/System/Library/Frameworks/OpenGL.framework/OpenGL" CACHE FILEPATH "opengl lib file")
  else()
    set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include" CACHE PATH "the include path for TBB")
    set(TBB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libtbb.so" "/usr/lib/x86_64-linux-gnu/libtbb.so.2" CACHE FILEPATH "the tbb .so files")
    set(OPENGL_LIBRARIES "/usr/lib/x86_64-linux-gnu/libOpenGL.so.0.0.0" CACHE FILEPATH "opengl lib file")
  endif(APPLE)
endif (UNIX)

set(STB
../stb/stb_image.h
../stb/stb_image_write.h
)

set(ICP
../icp/icp.h
../icp/icp.cpp
../icp/icp_point_to_plane.h
../icp/icp_point_to_plane.cpp
../icp/icp_point_to_point.h
../icp/icp_point_to_point.cpp
)

set(HDRS
camera.h
canvas.h
db.h
deform_tool.h
distance_map.h
ear_detector.h
face_detector.h
fill_holes.h
hdf5aux.h
im.h
io.h
jet.h
lscm.h
matcap.h
mesh.h
mm.h
morphable_model.h
mouse.h
objects.h
pc.h
pixel.h
primitives.h
scene.h
shape_predictor.h
sp.h
trackball.h
view.h
)
	
set(SRCS
meshscript.rc
camera.cpp
canvas.cpp
db.cpp
deform_tool.cpp
distance_map.cpp
ear_detector.cpp
face_detector.cpp
fill_holes.cpp
im.cpp
io.cpp
jet.cpp
lscm.cpp
main.cpp
matcap.cpp
mesh.cpp
mm.cpp
pc.cpp
pixel.cpp
scene.cpp
shape_predictor.cpp
sp.cpp
trackball.c
view.cpp
)

set(RESOURCES
meshscript.ico
)

if (WIN32)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:16000000")
endif (WIN32)

if (UNIX)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -msse4.1 -pthread")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse4.1 -pthread -std=c++11")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -msse4.1 -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -msse4.1 -pthread -std=c++11")
endif (UNIX)

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-DMEMORY_LEAK_TRACKING)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

include ("../jtk/jtk/jtk.cmake")

add_executable(meshscript ${HDRS} ${SRCS} ${RESOURCES} ${STB} ${ICP})

source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})
source_group("stb" FILES ${STB})
source_group("icp" FILES ${ICP})

target_include_directories(meshscript
    PRIVATE
    .
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../jtk/
    ${CMAKE_CURRENT_SOURCE_DIR}/../skiwi/
    ${CMAKE_CURRENT_SOURCE_DIR}/../stb/
    ${CMAKE_CURRENT_BINARY_DIR}/../hdf5
    ${CMAKE_CURRENT_BINARY_DIR}/../sdl2/include
    )	
	
target_link_libraries(meshscript
    PRIVATE
    ${OPENGL_LIBRARIES}    
    asm
    dlib
    libcork
    libpoisson
    libskiwi	
    trico
    hdf5-shared
    SDL2
    SDL2main
    )
    

if (UNIX)

set(XLIBINCLUDE "")
set(XLIBLIBRARY "X11.so")

if (APPLE)
set(XLIBINCLUDE "/usr/X11/include/")
set(XLIBLIBRARY "/usr/X11/lib/libX11.dylib")
endif (APPLE)

target_include_directories(meshscript
   PRIVATE
   ${XLIBINCLUDE}
   )
   
target_link_libraries(meshscript
    PRIVATE
    ${XLIBLIBRARY}
    )
endif (UNIX)

if (${JTK_THREADING} STREQUAL "tbb")
  target_include_directories(meshscript
      PRIVATE
      ${TBB_INCLUDE_DIR}
      )
  target_link_libraries(meshscript
      PRIVATE	
      ${TBB_LIBRARIES}      
      )	
endif (${JTK_THREADING} STREQUAL "tbb")

if (WIN32)
add_custom_command(TARGET meshscript POST_BUILD 
   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/../skiwi/bin" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif (WIN32)

add_custom_command(TARGET meshscript POST_BUILD 
   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/resources" "$<TARGET_FILE_DIR:meshscript>/resources")
   
add_custom_command(TARGET meshscript POST_BUILD 
   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/../skiwi/libskiwi/scm" "$<TARGET_FILE_DIR:meshscript>/scm")   
   
   
