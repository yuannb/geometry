cmake_minimum_required(VERSION 3.14)

project(GWB VERSION 1.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-fPIC)

SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

#将编译的汇编保存
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps=obj")

set(gwb_global_dir "${PROJECT_SOURCE_DIR}")
string(REGEX REPLACE "(.+)\\geometry.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})
set(pre_gwb_global_dir "${PROJECT_INIT_PATH}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Topo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Geo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tools)
add_subdirectory(${PROJECT_INIT_PATH}/3rd/googletest ${PROJECT_BINARY_DIR}/googleTest)

enable_testing()

list(APPEND EXTRA_LIBS topo)
list(APPEND EXTRA_LIBS geo)
list(APPEND EXTRA_LIBS tools)
list(APPEND EXTRA_LIBS gtest_main)

list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Topo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Geo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/tools")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/NurbsTools")
list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/eigen/include/eigen3")
list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")

add_executable(GWB nurbs_test.cpp debug_used.cpp)
add_executable(GWB1 test2.cpp)

target_link_libraries(GWB PUBLIC ${EXTRA_LIBS})
target_link_libraries(GWB1 PUBLIC ${EXTRA_LIBS})

target_include_directories(GWB PUBLIC  ${EXTRA_INCLUDES})
target_include_directories(GWB1 PUBLIC ${EXTRA_INCLUDES})

add_test(NAME nurbs_test COMMAND GWB)
