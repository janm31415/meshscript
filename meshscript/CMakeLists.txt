if (WIN32)
  set(SDL_INCLUDE_DIR "C:/Program Files/SDL2/include/SDL2/")
  set(SDL_LIBRARIES "C:/Program Files/SDL2/lib/SDL2.lib" "C:/Program Files/SDL2/lib/SDL2main.lib")
  set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include")
  set(TBB_LIBRARIES "C:/Program Files/TBB/lib/intel64/vc14/tbb.lib" "C:/Program Files/TBB/lib/intel64/vc14/tbb_debug.lib")
  set(OPENGL_LIBRARIES opengl32.lib)
endif (WIN32)

if (UNIX)
  set(SDL_INCLUDE_DIR "/usr/include/SDL2/")
  set(SDL_LIBRARIES "/usr/lib/x86_64-linux-gnu/libSDL2-2.0.so" "/usr/lib/x86_64-linux-gnu/libSDL2-2.0.so.0")
  set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include")
  set(TBB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libtbb.so" "/usr/lib/x86_64-linux-gnu/libtbb.so.2")
  set(OPENGL_LIBRARIES "/usr/lib/x86_64-linux-gnu/libOpenGL.so.0.0.0")
endif (UNIX)

if (APPLE)
  set(SDL_INCLUDE_DIR "/Library/Frameworks/SDL2.framework/Headers/")
  set(SDL_LIBRARIES "/Library/Frameworks/SDL2.framework/SDL2")
  set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include")
  set(TBB_LIBRARIES "/usr/local/lib/libtbb.dylib")
  set(OPENGL_LIBRARIES "/System/Library/Frameworks/OpenGL.framework/OpenGL")
endif (APPLE)

set(STB
../stb/stb_image.h
../stb/stb_image_write.h
)

set(PLY
#../ply/ply.h
#../ply/plyfile.c
../ply/rply.c
../ply/rply.h
../ply/rplyfile.h
)

set(ICP
../icp/icp.h
../icp/icp.cpp
../icp/icp_point_to_plane.h
../icp/icp_point_to_plane.cpp
../icp/icp_point_to_point.h
../icp/icp_point_to_point.cpp
../icp/point_tree.h
)

set(HDRS
camera.h
canvas.h
db.h
distance_map.h
ear_detector.h
face_detector.h
hashed_heap.h
hdf5aux.h
io.h
jet.h
matcap.h
mesh.h
mm.h
morphable_model.h
mouse.h
objects.h
pc.h
pixel.h
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
distance_map.cpp
ear_detector.cpp
face_detector.cpp
io.cpp
jet.cpp
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

add_executable(meshscript ${HDRS} ${SRCS} ${RESOURCES} ${STB} ${PLY} ${ICP})

source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})
source_group("stb" FILES ${STB})
source_group("ply" FILES ${PLY})
source_group("icp" FILES ${ICP})

target_include_directories(meshscript
    PRIVATE
    .
    ${CMAKE_CURRENT_SOURCE_DIR}/../
    ${CMAKE_CURRENT_SOURCE_DIR}/../jtk/
    ${CMAKE_CURRENT_SOURCE_DIR}/../skiwi/
    ${CMAKE_CURRENT_SOURCE_DIR}/../stb/
    ${CMAKE_CURRENT_SOURCE_DIR}/../ply/
    ${SDL_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}/../hdf5
    )	
	
target_link_libraries(meshscript
    PRIVATE
    ${OPENGL_LIBRARIES}
    ${SDL_LIBRARIES}
    asm
    dlib
    libpoisson
    libskiwi	
    hdf5-shared
    )	

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

#add_custom_command(TARGET meshscript POST_BUILD 
#   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/data" "$<TARGET_FILE_DIR:meshscript>/data")
   
add_custom_command(TARGET meshscript POST_BUILD 
   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/../skiwi/libskiwi/scm" "$<TARGET_FILE_DIR:meshscript>/scm")   
   
   
