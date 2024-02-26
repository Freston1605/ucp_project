#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "openbabel" for configuration "Release"
set_property(TARGET openbabel APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(openbabel PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libopenbabel.so.7.0.0"
  IMPORTED_SONAME_RELEASE "libopenbabel.so.7"
  )

list(APPEND _cmake_import_check_targets openbabel )
list(APPEND _cmake_import_check_files_for_openbabel "${_IMPORT_PREFIX}/lib/libopenbabel.so.7.0.0" )

# Import target "inchi" for configuration "Release"
set_property(TARGET inchi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(inchi PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libinchi.so.0.4.1"
  IMPORTED_SONAME_RELEASE "libinchi.so.0"
  )

list(APPEND _cmake_import_check_targets inchi )
list(APPEND _cmake_import_check_files_for_inchi "${_IMPORT_PREFIX}/lib/libinchi.so.0.4.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
