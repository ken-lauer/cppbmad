cmake_minimum_required(VERSION 3.15...3.27)

set(LIBNAME cppbmad)
set(PYBIND_LIBNAME _pybmad)

string(TOUPPER "${CMAKE_BUILD_TYPE}" _uppercase_cmake_build_type)
if(_uppercase_cmake_build_type STREQUAL "DEBUG")
  set(_default_acc_build_type "debug")
else()
  set(_default_acc_build_type "production")
endif()

set(BUILD_TYPE "${_default_acc_build_type}"
  CACHE STRING "Build type for ACC_ROOT_DIR paths (e.g., production or debug)"
)
set(CMAKE_EXPORT_COMPILE_COMMANDS true)

project(
  ${LIBNAME}
  VERSION 0.0.1
  LANGUAGES CXX Fortran
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(${LIBNAME} SHARED
  src/generated/proxy_mod.f90
  src/generated/bmad_routines.f90
  src/generated/sim_utils_routines.f90
  src/generated/tao_routines.f90
  src/tao_proxy_helpers.f90
  src/test_struct_defs.f90
  src/tests/routines/cppbmad_tests.f90
  src/generated/cppbmad_test_routines.f90

  src/common_structs.cpp
  src/generated/proxy.cpp
  src/generated/cppbmad_test_routines.cpp
  src/generated/bmad_routines.cpp
  src/generated/sim_utils_routines.cpp
  src/generated/tao_routines.cpp
  src/generated/to_string.cpp
  src/tao_index.cpp
  src/proxy_base.cpp
  src/to_string.cpp
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

if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  # No address sanitization options
else()
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address,undefined")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")
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

if(SKBUILD)
  install(TARGETS ${LIBNAME} DESTINATION pybmad)
else()
  install(TARGETS ${LIBNAME} DESTINATION lib)
endif()

foreach(TEST_NAME test_integration test_csr test_arrays)
  add_executable(${TEST_NAME} src/tests/${TEST_NAME}.cpp)
  target_include_directories(${TEST_NAME} PRIVATE "./include")
  target_link_directories(
    ${TEST_NAME}
    PRIVATE
    "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib"
  )
  target_link_libraries(${TEST_NAME} PRIVATE ${LIBNAME})
  set_property(TARGET ${TEST_NAME} PROPERTY CXX_STANDARD 17)
endforeach()

# Unit tests with doctest
add_executable(test_all_encompassing
  src/tests/test_main.cpp
  src/tests/test_all_encompassing_proxy.cpp
)
target_include_directories(test_all_encompassing PRIVATE
  "./include"
  "./src/tests"  # for doctest.h
)
target_link_directories(
  test_all_encompassing
  PRIVATE
  "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib"
)
target_link_libraries(test_all_encompassing PRIVATE ${LIBNAME})
set_property(TARGET test_all_encompassing PROPERTY CXX_STANDARD 17)

# Enable CTest
enable_testing()
add_test(NAME AllEncompassingProxyTests COMMAND test_all_encompassing)

option(BUILD_PYBMAD "Build the Bmad Python bindings" ON)

if(BUILD_PYBMAD OR SKBUILD)
  add_subdirectory(python)
endif()
