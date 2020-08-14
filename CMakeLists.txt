cmake_minimum_required(VERSION 3.12)
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

add_subdirectory(hdf5)
add_subdirectory(jtk)
add_subdirectory(meshscript)
add_subdirectory(skiwi)



set_property(GLOBAL PROPERTY USE_FOLDERS ON)


set_target_properties (asm PROPERTIES FOLDER skiwi)
set_target_properties (asm.tests PROPERTIES FOLDER skiwi)
set_target_properties (libskiwi PROPERTIES FOLDER skiwi)
set_target_properties (libskiwi.tests PROPERTIES FOLDER skiwi)
set_target_properties (life PROPERTIES FOLDER skiwi)
set_target_properties (s PROPERTIES FOLDER skiwi)
set_target_properties (utils PROPERTIES FOLDER skiwi)

set_target_properties (jtk PROPERTIES FOLDER jtk)
set_target_properties (jtk.tests PROPERTIES FOLDER jtk)

set_target_properties (hdf5-shared PROPERTIES FOLDER hdf5)
set_target_properties (hdf5-static PROPERTIES FOLDER hdf5)
set_target_properties (gen_hdf5-shared PROPERTIES FOLDER hdf5)
set_target_properties (gen_hdf5-static PROPERTIES FOLDER hdf5)
set_target_properties (H5make_libsettings PROPERTIES FOLDER hdf5)