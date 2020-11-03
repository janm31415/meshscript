
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
