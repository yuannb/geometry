
# include_directories(${gwb_global_dir}/include/tools)
# include_directories(${gwb_global_dir}/include/Topo)
# include_directories(${gwb_global_dir}/include/Geo)

# include_directories(${gwb_global_dir}/include/Geo)
file(GLOB INC_FILES ${gwb_global_dir}/include/tools/*.h)
file(GLOB SRC_FILES *.cpp)

# aux_source_directory(. SOURCES)

add_library(tools ${INC_FILES} ${SRC_FILES})
target_include_directories(tools PUBLIC ${gwb_global_dir}/include/tools)
target_include_directories(tools PUBLIC ${gwb_global_dir}/include/Topo)
target_include_directories(tools PUBLIC ${gwb_global_dir}/include/Geo)

target_link_libraries(tools Eigen3::Eigen)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
# source_group_by_dir_prefixed(src SRC_FILES)

# aux_source_directory(. SOURCES)

# add_library(tools ${SOURCES})
