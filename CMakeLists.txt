cmake_minimum_required(VERSION 3.10)
project (meshscript)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_PDB_OUTPUT_DIRECTORY     "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HDF5_EXTERNALLY_CONFIGURED 1)
set(HDF5_INSTALL_BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(HDF5_INSTALL_LIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(BUILD_TESTING OFF)

set(TRICO_SHARED "no")

add_subdirectory(dlib)
add_subdirectory(hdf5)

if (WIN32)
set(CMAKE_C_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_CXX_FLAGS_DEBUG "/W4 /MP /GF /RTCu /Od /MDd /Zi")
set(CMAKE_C_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "/W4 /MP /GF /O2 /Ob2 /Oi /Ot /MD /Zi /DNDEBUG")
endif (WIN32)

add_subdirectory(jtk)
add_subdirectory(libcork)
add_subdirectory(libpoisson)
add_subdirectory(meshscript)
add_subdirectory(skiwi)
add_subdirectory(schemerlicht)
add_subdirectory(trico/lz4)
add_subdirectory(trico/trico)
add_subdirectory(SDL2)


set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set_target_properties (libcork PROPERTIES FOLDER cork)
set_target_properties (libpoisson PROPERTIES FOLDER poisson)

set_target_properties (schemerlicht PROPERTIES FOLDER schemerlicht)
set_target_properties (schemerlicht.tests PROPERTIES FOLDER schemerlicht)
set_target_properties (sl PROPERTIES FOLDER schemerlicht)
set_target_properties (primdispatch PROPERTIES FOLDER schemerlicht/tools)
set_target_properties (primextract PROPERTIES FOLDER schemerlicht/tools)
set_target_properties (gameoflife PROPERTIES FOLDER schemerlicht/examples)


set_target_properties (asm PROPERTIES FOLDER skiwi)
set_target_properties (asm.tests PROPERTIES FOLDER skiwi)
set_target_properties (libskiwi PROPERTIES FOLDER skiwi)
set_target_properties (libskiwi.tests PROPERTIES FOLDER skiwi)
set_target_properties (life PROPERTIES FOLDER skiwi)
set_target_properties (s PROPERTIES FOLDER skiwi)

set_target_properties (dlib PROPERTIES FOLDER dlib)

set_target_properties (jtk PROPERTIES FOLDER jtk)
set_target_properties (jtk.tests PROPERTIES FOLDER jtk)
set_target_properties (jtk.static.tests PROPERTIES FOLDER jtk)

set_target_properties (trico PROPERTIES FOLDER trico)
set_target_properties (lz4 PROPERTIES FOLDER trico)

set_target_properties (hdf5-shared PROPERTIES FOLDER hdf5)
set_target_properties (hdf5-static PROPERTIES FOLDER hdf5)
set_target_properties (gen_hdf5-shared PROPERTIES FOLDER hdf5)
set_target_properties (gen_hdf5-static PROPERTIES FOLDER hdf5)
set_target_properties (H5make_libsettings PROPERTIES FOLDER hdf5)

set_target_properties (SDL2 PROPERTIES FOLDER SDL2)
set_target_properties (SDL2main PROPERTIES FOLDER SDL2)
set_target_properties (uninstall PROPERTIES FOLDER SDL2)
set_target_properties (sdl_headers_copy PROPERTIES FOLDER SDL2)
set_target_properties (SDL2_test PROPERTIES FOLDER SDL2)

