cmake_minimum_required(VERSION 3.14)
project(GWB VERSION 1.0)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# windows 编译
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")

add_definitions(-D_USE_MATH_DEFINES)

# add_compile_options(-fPIC)

#SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
#SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

#将编译的汇编保存
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -save-temps=obj")

set(gwb_global_dir "${PROJECT_SOURCE_DIR}")
string(REGEX REPLACE "(.+)\\geometry.*" "\\1" PROJECT_INIT_PATH ${PROJECT_SOURCE_DIR})
message(${PROJECT_INIT_PATH})
set(pre_gwb_global_dir "${PROJECT_INIT_PATH}")

add_subdirectory(${PROJECT_SOURCE_DIR}/src/Topo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/Geo)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tools)
add_subdirectory(${PROJECT_INIT_PATH}/3rd/googletest ${PROJECT_BINARY_DIR}/googleTest)
add_subdirectory(${PROJECT_INIT_PATH}/3rd/OpenMesh ${PROJECT_BINARY_DIR}/OpenMesh)
find_package(OpenMesh)

enable_testing()

list(APPEND EXTRA_LIBS topo)
list(APPEND EXTRA_LIBS geo)
list(APPEND EXTRA_LIBS tools)
list(APPEND EXTRA_LIBS gtest_main)

list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/eigen/include/eigen3")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Topo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/Geo")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/tools")
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/include/NurbsTools")
list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/libhv/cpputil")
list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/googletest/googletest/include")
list(APPEND EXTRA_INCLUDES "${PROJECT_INIT_PATH}/3rd/OpenMesh/src")

add_executable(GWB nurbs_test.cpp debug_used.cpp)
add_executable(ICURVET intersection_curve_test.cpp debug_used.cpp)
add_executable(KDTreeTest kdtree_test.cpp debug_used.cpp)
# big obj
if(CMAKE_COMPILER_IS_GNUCXX)
  # 检查GNU编译是否支持 -Wa,-mbig-obj 
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-Wa,-mbig-obj" GNU_BIG_OBJ_FLAG_ENABLE)
  message(STATUS GNU_BIG_OBJ_FLAG_ENABLE=${GNU_BIG_OBJ_FLAG_ENABLE})
endif()
target_compile_options(GWB 
	PRIVATE 
	$<$<CXX_COMPILER_ID:MSVC>:/bigobj> 
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<BOOL:${GNU_BIG_OBJ_FLAG_ENABLE}>>:-Wa,-mbig-obj>)
target_compile_options(ICURVET 
	PRIVATE 
	$<$<CXX_COMPILER_ID:MSVC>:/bigobj> 
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<BOOL:${GNU_BIG_OBJ_FLAG_ENABLE}>>:-Wa,-mbig-obj>)
target_compile_options(KDTreeTest 
	PRIVATE 
	$<$<CXX_COMPILER_ID:MSVC>:/bigobj> 
	$<$<AND:$<CXX_COMPILER_ID:GNU>,$<BOOL:${GNU_BIG_OBJ_FLAG_ENABLE}>>:-Wa,-mbig-obj>)

add_executable(GWB1 test2.cpp)
add_executable (OpenMeshTest build_cube_test.cpp)
target_link_libraries(OpenMeshTest PRIVATE OpenMeshCore OpenMeshTools)

target_link_libraries(GWB PUBLIC ${EXTRA_LIBS})
target_link_libraries(KDTreeTest PUBLIC ${EXTRA_LIBS})
target_link_libraries(GWB1 PUBLIC ${EXTRA_LIBS})
target_link_libraries(ICURVET PUBLIC ${EXTRA_LIBS})

target_include_directories(ICURVET PUBLIC  ${EXTRA_INCLUDES})
target_include_directories(KDTreeTest PUBLIC  ${EXTRA_INCLUDES})
target_include_directories(GWB1 PUBLIC ${EXTRA_INCLUDES})
target_include_directories(OpenMeshTest PUBLIC ${EXTRA_INCLUDES})
target_include_directories(ICURVET PUBLIC  ${EXTRA_INCLUDES})

add_test(NAME nurbs_test COMMAND GWB)
add_test(NAME ICurve_test COMMAND ICURVET)
add_test(NAME kdtree_test COMMAND KDTreeTest)
