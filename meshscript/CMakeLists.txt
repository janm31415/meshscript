
set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include")
set(TBB_LIBRARY "C:/Program Files/TBB/lib/intel64/vc14")
set(TBB_DIR "C:/Program Files/TBB")

find_package(TBB REQUIRED)

set(SDL_INCLUDE_DIR "C:/Program Files/SDL2/include")
set(SDL_LIBRARY "C:/Program Files/SDL2/lib")
find_package(SDL REQUIRED)

include_directories(${SDL_INCLUDE_DIR}/sdl2)
link_directories(${SDL_LIBRARY})

include_directories(${TBB_INCLUDE_DIR})
link_directories(${TBB_LIBRARY})

set(STB
../stb/stb_image.h
../stb/stb_image_write.h
)

set(HDRS
camera.h
canvas.h
db.h
jet.h
matcap.h
mesh.h
mouse.h
pixel.h
scene.h
trackball.h
view.h
)
	
set(SRCS
meshscript.rc
camera.cpp
canvas.cpp
db.cpp
jet.cpp
main.cpp
matcap.cpp
mesh.cpp
pixel.cpp
scene.cpp
trackball.c
view.cpp
)

set(RESOURCES
meshscript.ico
)

set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")

# general build definitions
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-DMEMORY_LEAK_TRACKING)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

add_executable(meshscript ${HDRS} ${SRCS} ${RESOURCES} ${STB})
source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})
source_group("stb" FILES ${STB})

target_include_directories(meshscript
    PRIVATE
    .
    ${CMAKE_CURRENT_SOURCE_DIR}/../jtk/
    ${CMAKE_CURRENT_SOURCE_DIR}/../skiwi/
    ${CMAKE_CURRENT_SOURCE_DIR}/../stb/
    )	
	
target_link_libraries(meshscript
    PRIVATE
    tbb.lib
    SDL2.lib
    SDL2main.lib
    opengl32.lib
    asm
    libskiwi	
    )	

if (WIN32)
add_custom_command(TARGET meshscript POST_BUILD 
   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/../skiwi/bin" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif (WIN32)