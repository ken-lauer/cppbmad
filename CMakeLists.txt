cmake_minimum_required(VERSION 3.15...3.27)

set(LIBNAME cppbmad)
set(BUILD_TYPE "debug" CACHE STRING "Build type for ACC_ROOT_DIR paths (e.g., production or debug)")

project(
  ${LIBNAME}
  VERSION 0.0.1
  LANGUAGES CXX Fortran
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB GENERATED_SOURCES src/generated/*.cpp)

add_library(${LIBNAME} SHARED
  ${GENERATED_SOURCES}
  src/_routines.cpp
  src/_routines.f90
  src/tao_index_proxies.cpp
  src/cpptest.cpp
)

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
  bmad
  sim_utils
  cpp_bmad_interface
)

if(APPLE)
  target_link_options(${LIBNAME} PRIVATE "-Wl,-no_compact_unwind")
endif()

set_property(TARGET ${LIBNAME} PROPERTY CXX_STANDARD 17)
target_compile_definitions(${LIBNAME} PRIVATE VERSION_INFO="${PROJECT_VERSION}")

install(TARGETS ${LIBNAME} DESTINATION lib)
