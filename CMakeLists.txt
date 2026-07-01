cmake_minimum_required(VERSION 3.16...4.3)
include(CMakePackageConfigHelpers)

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

if(DEFINED ENV{PKG_VERSION})
  set(PROJECT_VERSION_STRING "$ENV{PKG_VERSION}")
else()
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE PROJECT_VERSION_STRING
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    string(REPLACE "-" "." PROJECT_VERSION_STRING "${PROJECT_VERSION_STRING}")
  endif()
endif()
if(NOT PROJECT_VERSION_STRING)
  set(PROJECT_VERSION_STRING "0.0.0")
endif()

project(
  ${LIBNAME}
  VERSION ${PROJECT_VERSION_STRING}
  LANGUAGES CXX Fortran
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_library(${LIBNAME} SHARED
  src/generated/proxy_mod.f90
  src/generated/bmad_routines.f90
  src/generated/sim_utils_routines.f90
  src/generated/tao_routines.f90
  src/generated/bsim_routines.f90
  src/tao_proxy_helpers.f90
  src/test_struct_defs.f90
  src/tests/routines/cppbmad_tests.f90
  src/generated/cppbmad_test_routines.f90
  src/generated/cppbmad_extra_routines.f90
  src/array_desc.f90
  src/extra/ele_misalignments.f90

  src/common_structs.cpp
)

set(CPPBMAD_DOLLAR_OK_SOURCES
  src/extra/ele_misalignments.f90
)

# Re-open the library target so any source list below stays consistent.
target_sources(${LIBNAME} PRIVATE
  src/generated/proxy.cpp
  src/generated/cppbmad_test_routines.cpp
  src/generated/cppbmad_extra_routines.cpp
  src/generated/bmad_routines.cpp
  src/generated/sim_utils_routines.cpp
  src/generated/bsim_routines.cpp
  src/generated/tao_routines.cpp
  src/generated/to_string.cpp
  src/proxy_base.cpp
  src/to_string.cpp
)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(STATUS "Configuring gfortran options for ${LIBNAME}")

  set_source_files_properties(${CPPBMAD_DOLLAR_OK_SOURCES}
    PROPERTIES COMPILE_OPTIONS "-fdollar-ok"
  )

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

# if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
#   # No address sanitization options
# else()
#   set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fsanitize=address,undefined")
#   set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address,undefined")
#   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")
# endif()

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
  bsim
)

if(APPLE)
  target_link_options(${LIBNAME} PRIVATE "-Wl,-no_compact_unwind")
endif()


set_property(TARGET ${LIBNAME} PROPERTY CXX_STANDARD 17)
target_compile_definitions(${LIBNAME} PRIVATE VERSION_INFO="${PROJECT_VERSION}")

# Large reused header(s) should be precompiled:
target_precompile_headers(${LIBNAME} PRIVATE
  "$<$<COMPILE_LANGUAGE:CXX>:${CMAKE_CURRENT_SOURCE_DIR}/include/bmad/generated/proxy.hpp>"
)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

if(SKBUILD)
  install(TARGETS ${LIBNAME} DESTINATION pybmad)
else()
  install(TARGETS ${LIBNAME}
    EXPORT cppbmadTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
  )
  install(DIRECTORY include/
    DESTINATION include
  )
  write_basic_package_version_file(
    "cppbmadConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )
  # Library location info
  install(EXPORT cppbmadTargets
    FILE cppbmadTargets.cmake
    NAMESPACE cppbmad::
    DESTINATION lib/cmake/cppbmad
  )
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/cppbmadConfig.cmake"
    "include(\"\${CMAKE_CURRENT_LIST_DIR}/cppbmadTargets.cmake\")\n"
  )
  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/cppbmadConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cppbmadConfigVersion.cmake"
    DESTINATION lib/cmake/cppbmad
  )
endif()

# Unit tests with doctest
add_executable(cppbmad_tests
  src/tests/test_main.cpp
  src/tests/test_all_encompassing_proxy.cpp
  src/tests/test_arrays.cpp
  src/tests/test_csr.cpp
  src/tests/test_enums.cpp
  src/tests/test_to_string.cpp
  src/tests/test_exceptions.cpp
)
target_include_directories(cppbmad_tests PRIVATE
  "./include"
  "./src/tests"  # for doctest.h
)
target_link_directories(cppbmad_tests PRIVATE "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib")
target_link_libraries(cppbmad_tests PRIVATE ${LIBNAME})
set_property(TARGET cppbmad_tests PROPERTY CXX_STANDARD 17)

# Enable CTest
enable_testing()
add_test(NAME cppbmad_tests COMMAND cppbmad_tests)

# Python / pybmad
find_package(Python 3.10 COMPONENTS Interpreter Development.Module QUIET)

if(Python_FOUND)
  # nanobind stores its cmake config inside the Python package
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_VARIABLE _nanobind_cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _nanobind_result
  )
  if(_nanobind_result EQUAL 0 AND _nanobind_cmake_dir)
    find_package(nanobind QUIET PATHS "${_nanobind_cmake_dir}" NO_DEFAULT_PATH)
  endif()
endif()

if(nanobind_FOUND)
  set(_default_build_pybmad ON)
else()
  set(_default_build_pybmad OFF)
endif()
option(BUILD_PYBMAD "Build the Bmad Python bindings" ${_default_build_pybmad})

if(BUILD_PYBMAD OR SKBUILD)
  if(NOT Python_FOUND)
    find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  if(NOT nanobind_FOUND)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
      OUTPUT_VARIABLE _nanobind_cmake_dir
      OUTPUT_STRIP_TRAILING_WHITESPACE
      COMMAND_ERROR_IS_FATAL ANY
    )
    find_package(nanobind CONFIG REQUIRED PATHS "${_nanobind_cmake_dir}" NO_DEFAULT_PATH)
  endif()
  add_subdirectory(python)
endif()
