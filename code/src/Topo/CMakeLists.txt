
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb ")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")


file(GLOB INC_FILES ${gwb_global_dir}/include/Topo/*.h )
file(GLOB SRC_FILES *.cpp)

add_library(topo ${INC_FILES} ${SRC_FILES})
target_include_directories(topo PUBLIC ${gwb_global_dir}/include/Topo)
target_include_directories(topo PUBLIC ${gwb_global_dir}/include/Geo)
target_include_directories(topo PUBLIC ${gwb_global_dir}/include/NurbsTools)
target_include_directories(topo PUBLIC ${pre_gwb_global_dir}/3rd/eigen/include/eigen3)


source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})

