
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_compile_options(-fPIC)


file(GLOB INC_FILES ${ng_global_dir}/include/*.hpp)
file(GLOB SRC_FILES *.cpp)

message(${INC_FILES})

add_library(node_graph ${INC_FILES} ${SRC_FILES})

target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/3rd/eigen/include/eigen3)
target_include_directories(node_graph PUBLIC ${ng_global_dir}/include)
target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/geometry/include/NurbsTools)
target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/geometry/include/Geo)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
