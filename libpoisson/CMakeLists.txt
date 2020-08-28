find_package(OpenMP)

set(HDRS
libpoisson_api.h
poisson_reconstruction_screened.h
)


set(HDRS_POISSON_SCREENED
poisson_screened/Allocator.h
poisson_screened/Array.h
poisson_screened/BinaryNode.h
poisson_screened/BSplineData.h
poisson_screened/CmdLineParser.h
poisson_screened/Factor.h
poisson_screened/FunctionData.h
poisson_screened/Geometry.h
poisson_screened/Hash.h
poisson_screened/MarchingCubes.h
poisson_screened/MAT.h
poisson_screened/MemoryUsage.h
poisson_screened/MultiGridOctreeData.h
poisson_screened/Octree.h
poisson_screened/PlyVertexMini.h
poisson_screened/PointStream.h
poisson_screened/Polynomial.h
poisson_screened/PPolynomial.h
poisson_screened/SparseMatrix.h
poisson_screened/MyTime.h
poisson_screened/Vector.h
poisson_screened/Array.inl
poisson_screened/BSplineData.inl
poisson_screened/CmdLineParser.inl
poisson_screened/FunctionData.inl
poisson_screened/Geometry.inl
poisson_screened/MAT.inl
poisson_screened/MultiGridOctreeData.inl
poisson_screened/MultiGridOctreeData.IsoSurface.inl
poisson_screened/MultiGridOctreeData.SortedTreeNodes.inl
poisson_screened/Octree.inl
poisson_screened/PointStream.inl
poisson_screened/Polynomial.inl
poisson_screened/PPolynomial.inl
poisson_screened/SparseMatrix.inl
poisson_screened/Vector.inl
	)
	
set(SRCS
poisson_reconstruction_screened.cpp
)

set(SRCS_POISSON_SCREENED
poisson_screened/CmdLineParser.cpp
poisson_screened/Factor.cpp
poisson_screened/Geometry.cpp
poisson_screened/MarchingCubes.cpp
poisson_screened/PlyFile.cpp
#poisson_screened/Time.cpp
)

if (WIN32)
set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
endif (WIN32)

# general build definitions
add_definitions(-D_UNICODE)
add_definitions(-DUNICODE)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)
add_definitions(-D_SILENCE_STDEXT_HASH_DEPRECATION_WARNINGS)

add_library(libpoisson SHARED ${HDRS} ${SRCS}  ${HDRS_POISSON_SCREENED} ${SRCS_POISSON_SCREENED})

source_group("Header Files" FILES ${HDRS})
source_group("Source Files" FILES ${SRCS})
source_group("Poisson Screened\\Header Files" FILES ${HDRS_POISSON_SCREENED})
source_group("Poisson Screened\\Source Files" FILES ${SRCS_POISSON_SCREENED})

target_include_directories(libpoisson
    PRIVATE
  	${CMAKE_CURRENT_SOURCE_DIR}/../jtk/
    )	
	
target_link_libraries(libpoisson
    PRIVATE	
    )	

if(OpenMP_CXX_FOUND)
    target_link_libraries(libpoisson PUBLIC OpenMP::OpenMP_CXX)
    add_definitions(-D_OPENMP)
endif()