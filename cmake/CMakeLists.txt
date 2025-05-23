#[[
Copyright (C) 2025 Ryan L. Guy & Dennis Fassbaender

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
]]

cmake_minimum_required(VERSION 3.5)

function(set_output_directories DIR)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR} PARENT_SCOPE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${DIR} PARENT_SCOPE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${DIR} PARENT_SCOPE)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${DIR} PARENT_SCOPE)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${DIR} PARENT_SCOPE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${DIR} PARENT_SCOPE)
endfunction(set_output_directories)

set(BUILD_DEBUG OFF CACHE BOOL "Turn on to build a debug release")

# Disable distribute source (ON -> OFF) to speed up build. Remember to
# turn back on before distribution!
option(DISTRIBUTE_SOURCE "Include source code in addon" ON)

option(WITH_OPENCL "Compile project with OpenCL features" OFF)
option(WITH_ENGINE_TEST "Compile engine test executable" OFF)
option(WITH_PYFLUID_MODULE "Compile pyfluid Python module" OFF)
option(WITH_MIXBOX "Compile with Mixbox pigment mixing feature" OFF)
SET(LIBRARY_SUFFIX "" CACHE STRING "Optional suffix to add to the end of the library name")

if(APPLE)
    set(MACOSX_DEPLOYMENT_TARGET 10.10)
endif()

project(bl_flip_fluids)
set(CMAKE_BUILD_TYPE Release)

set(FLUIDENGINE_VERSION_TYPE_IS_STABLE_BUILD TRUE)
set(FLUIDENGINE_VERSION_MAJOR 1)
set(FLUIDENGINE_VERSION_MINOR 8)
set(FLUIDENGINE_VERSION_REVISION 3)
set(FLUIDENGINE_VERSION_DATE "23-MAR-2025")

if(FLUIDENGINE_VERSION_TYPE_IS_STABLE_BUILD)
    set(FLUIDENGINE_VERSION_TYPE_LABEL "GitHub Release")
    set(FLUIDENGINE_VERSION_TYPE_IS_STABLE_BUILD_PYTHON "True")
else()
    set(FLUIDENGINE_VERSION_TYPE_LABEL "Experimental")
    set(FLUIDENGINE_VERSION_TYPE_IS_STABLE_BUILD_PYTHON "False")
endif()

set(FLUIDENGINE_VERSION_LABEL "${FLUIDENGINE_VERSION_MAJOR}.${FLUIDENGINE_VERSION_MINOR}.${FLUIDENGINE_VERSION_REVISION} ${FLUIDENGINE_VERSION_TYPE_LABEL} ${FLUIDENGINE_VERSION_DATE}")
#if(BUILD_DEBUG)
#    set(FLUIDENGINE_VERSION_LABEL "${FLUIDENGINE_VERSION_LABEL} (DEBUG BUILD)") 
#endif()

message(STATUS "FLIP Fluids version ${FLUIDENGINE_VERSION_LABEL}")
if(BUILD_DEBUG)
    message(STATUS "Building in debug mode")
else()
    message(STATUS "Building in release mode")
endif()

# Configure Mixbox
if(WITH_MIXBOX)
    set(FLUIDENGINE_IS_MIXBOX_SUPPORTED "True")
    add_definitions(-DWITH_MIXBOX=1)
    set(MIXBOX_SOURCE_CPP "src/engine/mixbox/mixbox.cpp")
else()
    set(FLUIDENGINE_IS_MIXBOX_SUPPORTED "False")
    add_definitions(-DWITH_MIXBOX=0)
    set(MIXBOX_SOURCE_CPP "src/engine/mixbox/mixbox_stub.cpp")
endif()

# Configure OpenCL
if(WITH_OPENCL)
    find_package(OpenCL)
    if(NOT OpenCL_FOUND)
        message(WARNING "WARNING: OpenCL was not found on your system.\nInstall an OpenCL SDK specific to your GPU vender (AMD, NVIDIA, Intel, etc.). This project will be built with OpenCL features disabled.")
        set(WITH_OPENCL OFF)
    endif()
endif()

if(WITH_OPENCL)
    if(MSVC OR MINGW)
        # Set to OpenCL ICD loader instead of vendor library
        set(OpenCL_LIBRARY "C:/Windows/System32/OpenCL.dll")
    endif()

    # Allows use os clCreateCommandQueue in OpenCL versions >= 2.0
    add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_2_APIS)

    add_definitions(-DWITH_OPENCL=1)
    include_directories(src ${OpenCL_INCLUDE_DIRS})
else()
    add_definitions(-DWITH_OPENCL=0)
endif()

# Configure compiler/OS specific flags
if(MSVC)
    message(FATAL_ERROR "Error: Compilation using MSVC (Microsoft Visual Studio) is not supported. Building with MSVC will result in errors, performance issues, and broken simulation features. Change this FATAL_ERROR to WARNING in the CMake script to continue at your own risk.")
    add_definitions(-DCL_USE_DEPRECATED_OPENCL_1_1_APIS -DNOMINMAX -D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /W4 /EHsc /wd4244")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /LARGEADDRESSAWARE")
elseif(MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static -static-libgcc -static-libstdc++ -O3 -std=c++20 -Wall")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-stdcall-fixup")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--enable-stdcall-fixup")
elseif(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++20 -Wall -Wno-return-type-c-linkage")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_NANOSLEEP -fopenmp -O3 -std=c++20 -Wall -fPIC -pthread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${PROJECT_SOURCE_DIR}/src/engine/glibc_version_header/force_link_glibc_2.5.h")
endif()

if(BUILD_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_DEBUG")
endif()

# Configure Version
configure_file(
  "${PROJECT_SOURCE_DIR}/src/engine/versionutils.cpp.in"
  "${PROJECT_SOURCE_DIR}/src/engine/versionutils.cpp"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/addon/__init__.py.in"
  "${PROJECT_SOURCE_DIR}/src/addon/__init__.py"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/addon/blender_manifest.toml.in"
  "${PROJECT_SOURCE_DIR}/src/addon/blender_manifest.toml"
)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/addon/utils/installation_utils.py.in"
  "${PROJECT_SOURCE_DIR}/src/addon/utils/installation_utils.py"
)

# Generate kernels.cpp
file(READ "src/engine/kernels/trilinearinterpolate.cl" TRILINEAR_INTERPOLATE_KERNEL_STRING)
file(READ "src/engine/kernels/scalarfield.cl" SCALAR_FIELD_KERNEL_STRING)
configure_file(
  "${PROJECT_SOURCE_DIR}/src/engine/kernels/kernels.cpp.in"
  "${PROJECT_SOURCE_DIR}/src/engine/kernels/kernels.cpp"
)

# Object library
file(GLOB SOURCES "src/engine/*.cpp" "src/engine/c_bindings/*.cpp" "src/engine/kernels/*.cpp" "src/engine/opencl_bindings/*.cpp" ${MIXBOX_SOURCE_CPP})
add_library(objects OBJECT ${SOURCES})

if(WITH_ENGINE_TEST)
    # Test executable
    set(TEST_DIR "${CMAKE_BINARY_DIR}/bl_flip_fluids/test")
    set_output_directories("${TEST_DIR}")
    add_executable(bl_flip_fluids $<TARGET_OBJECTS:objects>)
    set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "engine_test")
    if(WITH_OPENCL)
        target_link_libraries(bl_flip_fluids "${OpenCL_LIBRARY}")
    endif()
endif()

if(WITH_PYFLUID_MODULE)
    # Pyfluid library
    set(PYTHON_MODULE_DIR "${CMAKE_BINARY_DIR}/bl_flip_fluids/pyfluid")
    set(PYTHON_MODULE_LIB_DIR "${CMAKE_BINARY_DIR}/bl_flip_fluids/pyfluid/lib")
    set_output_directories("${PYTHON_MODULE_LIB_DIR}")
    add_library(pyfluid SHARED $<TARGET_OBJECTS:objects>)
    if(WITH_OPENCL) 
        target_link_libraries(pyfluid "${OpenCL_LIBRARY}")
    endif()
    file(COPY "${CMAKE_SOURCE_DIR}/src/engine/pyfluid/" DESTINATION "${PYTHON_MODULE_DIR}")
endif()

# Blender addon
set(BLENDER_ADDON_DIR "${CMAKE_BINARY_DIR}/bl_flip_fluids/flip_fluids_addon")
set_output_directories("${BLENDER_ADDON_DIR}/pyfluid/lib")
if(BUILD_DEBUG)
    add_library("blpyfluiddebug${LIBRARY_SUFFIX}" SHARED $<TARGET_OBJECTS:objects>)
else()
    add_library("blpyfluidrelease${LIBRARY_SUFFIX}" SHARED $<TARGET_OBJECTS:objects>)
endif()

if(WITH_OPENCL) 
    if(BUILD_DEBUG)
        target_link_libraries("blpyfluiddebug${LIBRARY_SUFFIX}" "${OpenCL_LIBRARY}")
    else()
        target_link_libraries("blpyfluidrelease${LIBRARY_SUFFIX}" "${OpenCL_LIBRARY}")
    endif()
endif()

file(COPY "${CMAKE_SOURCE_DIR}/src/addon/" DESTINATION "${BLENDER_ADDON_DIR}")
file(COPY "${CMAKE_SOURCE_DIR}/src/engine/pyfluid/" DESTINATION "${BLENDER_ADDON_DIR}/pyfluid")

if(DISTRIBUTE_SOURCE)
    set(SOURCES_ARCHIVE_ROOT
        "${CMAKE_SOURCE_DIR}/__init__.py"
        "${CMAKE_SOURCE_DIR}/build.py"
        "${CMAKE_SOURCE_DIR}/.gitignore"
        "${CMAKE_SOURCE_DIR}/LICENSE_GPLv3.md"
        "${CMAKE_SOURCE_DIR}/LICENSE_MIT.md"
        "${CMAKE_SOURCE_DIR}/LICENSE_Standard_Royalty_Free.md"
        "${CMAKE_SOURCE_DIR}/README.md"
    )

    file(GLOB_RECURSE SOURCES_ARCHIVE "${CMAKE_SOURCE_DIR}/src/*" "${CMAKE_SOURCE_DIR}/cmake/*")
    list (APPEND SOURCES_ARCHIVE ${SOURCES_ARCHIVE_ROOT})
    file(MAKE_DIRECTORY ${BLENDER_ADDON_DIR}/source)
    file(ARCHIVE_CREATE OUTPUT "${BLENDER_ADDON_DIR}/source/flip_fluids_addon_source.zip" PATHS ${SOURCES_ARCHIVE} FORMAT "zip")
endif()
