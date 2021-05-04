set(KIM-API_MIN_VERSION 2.1.3)
find_package(CURL)
if(CURL_FOUND)
  if(CMAKE_VERSION VERSION_LESS 3.12)
    target_include_directories(lammps PRIVATE ${CURL_INCLUDE_DIRS})
    target_link_libraries(lammps PRIVATE ${CURL_LIBRARIES})
  else()
    target_link_libraries(lammps PRIVATE CURL::libcurl)
  endif()
  target_compile_definitions(lammps PRIVATE -DLMP_KIM_CURL)
  set(LMP_DEBUG_CURL OFF CACHE STRING "Set libcurl verbose mode on/off. If on, it displays a lot of verbose information about its operations.")
  mark_as_advanced(LMP_DEBUG_CURL)
  if(LMP_DEBUG_CURL)
    target_compile_definitions(lammps PRIVATE -DLMP_DEBUG_CURL)
  endif()
  set(LMP_NO_SSL_CHECK OFF CACHE STRING "Tell libcurl to not verify the peer. If on, the connection succeeds regardless of the names in the certificate. Insecure - Use with caution!")
  mark_as_advanced(LMP_NO_SSL_CHECK)
  if(LMP_NO_SSL_CHECK)
    target_compile_definitions(lammps PRIVATE -DLMP_NO_SSL_CHECK)
  endif()
endif()
find_package(PkgConfig QUIET)
set(DOWNLOAD_KIM_DEFAULT ON)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(KIM-API QUIET libkim-api>=${KIM-API_MIN_VERSION})
  if(KIM-API_FOUND)
    set(DOWNLOAD_KIM_DEFAULT OFF)
  endif()
endif()
option(DOWNLOAD_KIM "Download KIM-API from OpenKIM instead of using an already installed one" ${DOWNLOAD_KIM_DEFAULT})
if(DOWNLOAD_KIM)
  message(STATUS "KIM-API download requested - we will build our own")
  include(ExternalProject)
  enable_language(C)
  enable_language(Fortran)
  ExternalProject_Add(kim_build
    URL https://s3.openkim.org/kim-api/kim-api-2.1.3.txz
    URL_MD5 6ee829a1bbba5f8b9874c88c4c4ebff8
    BINARY_DIR build
    CMAKE_ARGS ${CMAKE_REQUEST_PIC}
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
               -DCMAKE_INSTALL_LIBDIR=lib
               -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
               -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
               BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libkim-api${CMAKE_SHARED_LIBRARY_SUFFIX}
    )
  ExternalProject_get_property(kim_build INSTALL_DIR)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include/kim-api)
  add_library(LAMMPS::KIM UNKNOWN IMPORTED)
  set_target_properties(LAMMPS::KIM PROPERTIES
    IMPORTED_LOCATION "${INSTALL_DIR}/lib/libkim-api${CMAKE_SHARED_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES "${INSTALL_DIR}/include/kim-api")
  target_link_libraries(lammps PRIVATE LAMMPS::KIM)
  add_dependencies(LAMMPS::KIM kim_build)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(KIM-API REQUIRED IMPORTED_TARGET libkim-api>=${KIM-API_MIN_VERSION})
  target_link_libraries(lammps PRIVATE PkgConfig::KIM-API)
endif()
