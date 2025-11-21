cmake_minimum_required(VERSION 3.15...3.27)

set(LIBNAME cppbmad)
set(BUILD_TYPE "debug" CACHE STRING "Build type for ACC_ROOT_DIR paths (e.g., production or debug)")

project(
  ${LIBNAME}
  VERSION 0.0.1
  LANGUAGES CXX Fortran
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(${LIBNAME} SHARED
  src/tao_proxy_helpers.f90
  src/generated/proxy_mod.f90
  src/test_struct_defs.f90

  src/proxy_base.cpp
  src/to_string.cpp
  src/generated/tao_proxies.cpp
)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Configuring gfortran options for ${LIBNAME}")

  # This strange syntax add compiler options just for Fortran source files:
  target_compile_options(${LIBNAME} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:
    # Generate backtrace on runtime error
    -fbacktrace

    # Allow lines to be any length (standard is usually 132)
    -ffree-line-length-none
  # Disable range checking on array subscripts
  # -fno-range-check

  # Allow dollar signs in entity names
  # -fdollar-ok

  # Generate backtrace on runtime error
  -fbacktrace

  # Link static libraries
  # add_link_options(-Bstatic

  # Allow lines to be any length (standard is usually 132
  -ffree-line-length-none

  # Enable OpenMP shared-memory parallel processing
  # -fopenmp
  # link options, too: -fopenmp
  >
  )

endif()

target_include_directories(${LIBNAME} PRIVATE
  "./include"
  "$ENV{ACC_ROOT_DIR}/cpp_bmad_interface/include"
  "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/modules"
)

target_link_directories(
  ${LIBNAME}
  PRIVATE
  "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib"
)

target_link_libraries(${LIBNAME} PUBLIC
  forest
  bmad
  sim_utils
  tao
  # cpp_bmad_interface
)

if(APPLE)
  target_link_options(${LIBNAME} PRIVATE "-Wl,-no_compact_unwind")
endif()

set_property(TARGET ${LIBNAME} PROPERTY CXX_STANDARD 17)
target_compile_definitions(${LIBNAME} PRIVATE VERSION_INFO="${PROJECT_VERSION}")

install(TARGETS ${LIBNAME} DESTINATION lib)

add_executable(test1 src/tests/test1.cpp)
target_include_directories(test1 PRIVATE
  "./include"
)
target_link_directories(
  test1
  PRIVATE
  "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib"
)
target_link_libraries(test1 PRIVATE ${LIBNAME})
