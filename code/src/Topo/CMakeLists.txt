
# include_directories(${gwb_global_dir}/include/Topo)


file(GLOB INC_FILES ${gwb_global_dir}/include/Topo/*.h )
file(GLOB SRC_FILES *.cpp)

add_library(topo ${INC_FILES} ${SRC_FILES})
target_include_directories(topo PUBLIC ${gwb_global_dir}/include/Topo)
target_link_libraries(topo Eigen3::Eigen)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
# source_group_by_dir_prefixed(src SRC_FILES)

# aux_source_directory(. SOURCES)

# add_library(topo ${SOURCES})
