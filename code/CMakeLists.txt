# TODO 1: Set the minimum required version of CMake to be 3.10
cmake_minimum_required(VERSION 3.10)

# set(CMAKE_CXX_COMPILER "C:/msys64/mingw64/bin/g++.exe")

# set(CMAKE_C_COMPILER "C:/msys64/mingw64/bin/gcc.exe")
# TODO 2: Create a project named Tutorial
project(GWB VERSION 1.0)

# TODO 7: Set the project version number as 1.0 in the above project command

# TODO 6: Set the variable CMAKE_CXX_STANDARD to 11
#         and the variable CMAKE_CXX_STANDARD_REQUIRED to True
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")


set(gwb_global_dir "${PROJECT_SOURCE_DIR}")

set(Eigen3_DIR C:/Users/Administrator/Desktop/workbench/3rd/eigen/share/eigen3/cmake)
find_package(Eigen3 REQUIRED NO_MODULE)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Topo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Geo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tools)
list(APPEND EXTRA_LIBS topo)
list(APPEND EXTRA_LIBS geo)
list(APPEND EXTRA_LIBS tools)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Topo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Geo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/tools")
list(APPEND EXTRA_INCLUDES "C:/Users/Administrator/Desktop/workbench/3rd/include")


aux_source_directory(. SOURCES)

# TODO 3: Add an executable called Tutorial to the project
# Hint: Be sure to specify the source file as tutorial.cxx
# add_executable(GWB params.cpp gwballocation.cpp query.cpp euler.cpp delaunay.cpp build.cpp )
add_executable(GWB ${SOURCES})

target_link_libraries(GWB PUBLIC ${EXTRA_LIBS})

target_include_directories(GWB PUBLIC "${PROJECT_SOURCE_DIR}"
                                        ${EXTRA_INCLUDES})
# TODO 9: Use target_include_directories to include ${PROJECT_BINARY_DIR}
# target_include_directories(GWB PUBLIC
#                             "${PROJECT_BINARY_DIR}")
