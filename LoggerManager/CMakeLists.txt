cmake_minimum_required(VERSION 3.14)

project(log VERSION 1.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-fPIC)

SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

set(ng_global_dir "${PROJECT_SOURCE_DIR}")
string(REGEX REPLACE "(.+)\\LoggerManager.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})
set(pre_ng_global_dir "${PROJECT_INIT_PATH}")

message("logerdir : ${ng_global_dir}/3rd/spdlog/include")
list(APPEND EXTRA_INCLUDES "${pre_ng_global_dir}/3rd/spdlog/include")

file(GLOB INC_FILES *.hpp)
file(GLOB SRC_FILES *.cpp)

add_library(log ${INC_FILES} ${SRC_FILES})
target_include_directories(log PUBLIC  ${EXTRA_INCLUDES})

if(MSVC)
    target_compile_options(log PRIVATE /utf-8)
endif()


source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
