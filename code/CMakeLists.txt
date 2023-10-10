# TODO 1: Set the minimum required version of CMake to be 3.10
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
set(third_package_dir "../3rd")

find_package(Eigen3 REQUIRED NO_MODULE)

# add_subdirectory(../3rd/thread-pool)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Topo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Geo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tools)
enable_testing()
string(REGEX REPLACE "(.+)\\code.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})

add_subdirectory(${PROJECT_INIT_PATH}/3rd/googletest ${PROJECT_BINARY_DIR}/googleTest)

# list(APPEND EXTRA_LIBS thread_pool)
list(APPEND EXTRA_LIBS topo)
list(APPEND EXTRA_LIBS geo)
list(APPEND EXTRA_LIBS tools)
list(APPEND EXTRA_LIBS gtest_main)
# list(APPEND EXTRA_INCLUDES "../3rd/thread_pool/include")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Topo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Geo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/tools")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/NurbsTools")


# list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")



aux_source_directory(. SOURCES)


add_executable(GWB nurbs_test.cpp)
add_executable(GWB1 test2.cpp)

target_link_libraries(GWB PUBLIC ${EXTRA_LIBS})
target_link_libraries(GWB1 PUBLIC ${EXTRA_LIBS})

target_include_directories(GWB PUBLIC "${PROJECT_SOURCE_DIR}"
                                        ${EXTRA_INCLUDES})

target_include_directories(GWB1 PUBLIC "${PROJECT_SOURCE_DIR}"
${EXTRA_INCLUDES})

target_include_directories(GWB PUBLIC "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")
target_include_directories(GWB1 PUBLIC "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")


add_test(NAME nurbs_test COMMAND GWB)



