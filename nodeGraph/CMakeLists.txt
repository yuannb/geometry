cmake_minimum_required(VERSION 3.14)

project(ng VERSION 1.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-fPIC)

SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

#将编译的汇编保存
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps=obj")

set(ng_global_dir "${PROJECT_SOURCE_DIR}")
string(REGEX REPLACE "(.+)\\nodeGraph.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})
set(pre_ng_global_dir "${PROJECT_INIT_PATH}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src)


# list(APPEND EXTRA_LIBS node_graph)
# list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include")


# list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/eigen/include/eigen3")
# list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/geometry/include/Geo")



# add_library(ng ReadJson.hpp ReadJson.cpp)

# target_link_libraries(ng PUBLIC ${EXTRA_LIBS})
# target_include_directories(ng PUBLIC  ${EXTRA_INCLUDES})