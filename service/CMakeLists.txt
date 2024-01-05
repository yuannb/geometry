cmake_minimum_required(VERSION 3.14)
project(service VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-fPIC)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

enable_testing()
string(REGEX REPLACE "(.+)\\service.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})
set(pre_global_dir "${PROJECT_INIT_PATH}")

add_subdirectory(${PROJECT_INIT_PATH}/3rd/googletest ${PROJECT_BINARY_DIR}/googleTest)
add_subdirectory(${PROJECT_INIT_PATH}/nodeGraph ${PROJECT_BINARY_DIR}/nodeGraph)
add_subdirectory(${PROJECT_INIT_PATH}/3rd/libhv ${PROJECT_BINARY_DIR}/libhv)

list(APPEND EXTRA_LIBS gtest_main)
list(APPEND EXTRA_LIBS node_graph)
list(APPEND EXTRA_LIBS hv_static)

list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include")
list(APPEND EXTRA_INCLUDES "${pre_global_dir}/nodeGraph/include")
list(APPEND EXTRA_INCLUDES "${pre_global_dir}/3rd/libhv/include/hv")

add_executable(service ReadJson.cpp test.cpp)

target_link_libraries(service PUBLIC ${EXTRA_LIBS})
target_include_directories(service PUBLIC "${EXTRA_INCLUDES}")

