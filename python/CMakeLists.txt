cmake_minimum_required(VERSION 3.16...4.3)

set(CMAKE_EXPORT_COMPILE_COMMANDS true)
set(PYBIND_LIBNAME _pybmad)

file(GLOB GENERATED_SOURCES
    CONFIGURE_DEPENDS
    "src/generated/*.cpp"
)

nanobind_add_module(${PYBIND_LIBNAME}
    ${GENERATED_SOURCES}
    src/arrays.cpp
    src/common_structs.cpp
)

# Include directories
target_include_directories(${PYBIND_LIBNAME} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    # cppbmad include:
    "${CMAKE_CURRENT_SOURCE_DIR}/../include"
    "$ENV{ACC_ROOT_DIR}/cpp_bmad_interface/include"
    "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/modules"
)

# Link directories
target_link_directories(${PYBIND_LIBNAME} PRIVATE "$ENV{ACC_ROOT_DIR}/${BUILD_TYPE}/lib")

# Link libraries
target_link_libraries(${PYBIND_LIBNAME} PRIVATE
    cppbmad
    bmad
    sim_utils
)

target_compile_definitions(${PYBIND_LIBNAME} PRIVATE VERSION_INFO=${PROJECT_VERSION})

# Large reused header(s) should be precompiled:
target_precompile_headers(${PYBIND_LIBNAME} PRIVATE
  <nanobind/nanobind.h>
  <nanobind/stl/array.h>
  <nanobind/stl/string.h>
  <nanobind/stl/vector.h>
  <nanobind/stl/optional.h>
  <nanobind/stl/complex.h>
  <nanobind/ndarray.h>
  "${CMAKE_CURRENT_SOURCE_DIR}/../include/bmad/generated/proxy.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/../include/bmad/generated/to_string.hpp"
)

# --- Installation Strategy for pyproject.toml ---

# When using scikit-build-core, we install into the python package directory
# relative to site-packages.
# logic: We want _pybmad.so (.dylib) to sit inside the `pybmad` folder alongside __init__.py.
install(TARGETS ${PYBIND_LIBNAME} DESTINATION pybmad)
# # Check if generated files exist using a marker file, or force run
# find_file(CODEGEN_MARKER "proxy.cpp" PATHS src/generated NO_DEFAULT_PATH)
#
# if(NOT CODEGEN_MARKER)
#   message(STATUS "Running python codegen step...")
#   execute_process(
#         COMMAND ${Python_EXECUTABLE} -m codegen
#         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#         RESULT_VARIABLE CODEGEN_RESULT
#     )
#   if(NOT CODEGEN_RESULT EQUAL 0)
#     message(FATAL_ERROR "Codegen failed")
#   endif()
# endif()

# Update RPATH so _pybmad.so finds libcppbmad.dylib/so in the same folder
if(APPLE)
  set_target_properties(${PYBIND_LIBNAME} PROPERTIES
    INSTALL_RPATH "@loader_path"
  )
else()
  set_target_properties(${PYBIND_LIBNAME} PROPERTIES
    INSTALL_RPATH "$ORIGIN"
  )
endif()

# Stub generation using nanobind's built-in stubgen
set(EXISTING_PYTHONPATH "$ENV{PYTHONPATH}")

if(EXISTING_PYTHONPATH)
  set(PYTHONPATH "$<TARGET_FILE_DIR:${PYBIND_LIBNAME}>:${EXISTING_PYTHONPATH}")
else()
  set(PYTHONPATH "$<TARGET_FILE_DIR:${PYBIND_LIBNAME}>")
endif()

add_custom_command(TARGET ${PYBIND_LIBNAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E env
  "PYTHONPATH=${PYTHONPATH}"
  "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:cppbmad>:$ENV{LD_LIBRARY_PATH}" # Linux
  "DYLD_LIBRARY_PATH=$<TARGET_FILE_DIR:cppbmad>:$ENV{DYLD_LIBRARY_PATH}" # MacOS

  "${Python_EXECUTABLE}" -m nanobind.stubgen
  -m "${PYBIND_LIBNAME}"
  -o "${CMAKE_CURRENT_SOURCE_DIR}/pybmad/${PYBIND_LIBNAME}.pyi"

  COMMENT "Generating Python stubs for ${PYBIND_LIBNAME}"
  VERBATIM
)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/pybmad/"
        DESTINATION pybmad
        FILES_MATCHING
        PATTERN "*.py"
        PATTERN "*.pyi"
        REGEX "__pycache__" EXCLUDE
        PATTERN "*.pyc" EXCLUDE
)
