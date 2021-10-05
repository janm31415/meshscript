if (WIN32)
  set(TBB_INCLUDE_DIR "C:/Program Files/TBB/include" CACHE PATH "the include path for TBB")
  set(TBB_LIBRARIES "C:/Program Files/TBB/lib/intel64/vc14/tbb.lib" "C:/Program Files/TBB/lib/intel64/vc14/tbb_debug.lib" CACHE FILEPATH "the tbb lib files")  
endif (WIN32)

if (UNIX)
  if (APPLE)
    set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include" CACHE PATH "the include path for TBB")
    set(TBB_LIBRARIES "/usr/local/lib/libtbb.dylib" CACHE FILEPATH "the tbb .dylib file")    
  else()
    set(TBB_INCLUDE_DIR "/usr/local/include" "/usr/include" CACHE PATH "the include path for TBB")
    set(TBB_LIBRARIES "/usr/lib/x86_64-linux-gnu/libtbb.so" "/usr/lib/x86_64-linux-gnu/libtbb.so.2" CACHE FILEPATH "the tbb .so files")    
  endif(APPLE)
endif (UNIX)

set(HDRS
aabvh.h
cork.h
union_find.h
fixint.h
fixisct.h
ga4.h
isct.h
libcork_api.h
quantization.h
triangle.h
)

set(SRCS
cork.cpp
triangle.c
)

# general build definitions
add_definitions(-D_UNICODE)
add_definitions(-DUNICODE)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
include ("../jtk/jtk/jtk.cmake")

add_library(libcork SHARED ${HDRS} ${SRCS})

source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})

target_include_directories(libcork
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../jtk/
    )	
	
target_link_libraries(libcork
    PRIVATE	
    )	


if (${JTK_THREADING} STREQUAL "tbb")
  target_include_directories(libcork
      PRIVATE
      ${TBB_INCLUDE_DIR}
      )
  target_link_libraries(libcork
      PRIVATE	
      ${TBB_LIBRARIES}      
      )	
endif (${JTK_THREADING} STREQUAL "tbb")