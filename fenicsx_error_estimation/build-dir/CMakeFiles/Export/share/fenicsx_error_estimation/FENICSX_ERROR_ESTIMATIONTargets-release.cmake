#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fenicsx_error_estimation" for configuration "Release"
set_property(TARGET fenicsx_error_estimation APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(fenicsx_error_estimation PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libfenicsx_error_estimation.so.0.5.0"
  IMPORTED_SONAME_RELEASE "libfenicsx_error_estimation.so.0.5"
  )

list(APPEND _IMPORT_CHECK_TARGETS fenicsx_error_estimation )
list(APPEND _IMPORT_CHECK_FILES_FOR_fenicsx_error_estimation "${_IMPORT_PREFIX}/lib/libfenicsx_error_estimation.so.0.5.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
