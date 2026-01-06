cmake_minimum_required(VERSION 3.15...3.27)

set(PYBIND_LIBNAME _pybmad)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

file(GLOB GENERATED_SOURCES
    CONFIGURE_DEPENDS
    "src/generated/*.cpp"
)

python_add_library(${PYBIND_LIBNAME} MODULE
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
# Note: "cppbmad" is visible here because it was defined in the parent CMake scope.
target_link_libraries(${PYBIND_LIBNAME} PUBLIC
    pybind11::headers
    cppbmad
    bmad
    sim_utils
    cpp_bmad_interface
)

target_compile_definitions(${PYBIND_LIBNAME} PRIVATE VERSION_INFO=${PROJECT_VERSION})

# --- Installation Strategy for pyproject.toml ---

# When using scikit-build-core, we install into the python package directory
# relative to site-packages.
# logic: We want _pybmad.so (.pyd) to sit inside the `pybmad` folder alongside __init__.py.
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
if (APPLE)
  set_target_properties(${PYBIND_LIBNAME} PROPERTIES
        INSTALL_RPATH "@loader_path"
    )
else()
  set_target_properties(${PYBIND_LIBNAME} PROPERTIES
        INSTALL_RPATH "$ORIGIN"
    )
endif()

# Only run this if we have the pybind11-stubgen tool available
find_program(STUBGEN pybind11-stubgen)

if(STUBGEN)
  add_custom_command(TARGET ${PYBIND_LIBNAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E env
      "PYTHONPATH=$<TARGET_FILE_DIR:${PYBIND_LIBNAME}>"
      "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:cppbmad>:$ENV{LD_LIBRARY_PATH}" # Linux
      "DYLD_LIBRARY_PATH=$<TARGET_FILE_DIR:cppbmad>:$ENV{DYLD_LIBRARY_PATH}" # MacOS
      ${STUBGEN} -o "${CMAKE_CURRENT_SOURCE_DIR}/pybmad" "${PYBIND_LIBNAME}"

      COMMENT "Generating Python Stubs for ${PYBIND_LIBNAME}"
      VERBATIM
  )

  # Install the generated stubs (.pyi files)
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/pybmad/"
            DESTINATION pybmad
            FILES_MATCHING PATTERN "*.pyi")
endif()
