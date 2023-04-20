
# include_directories(${gwb_global_dir}/include/Geo)
file(GLOB INC_FILES ${gwb_global_dir}/include/Geo/*.h)
file(GLOB SRC_FILES *.cpp)

# aux_source_directory(. SOURCES)

add_library(geo ${INC_FILES} ${SRC_FILES})
target_include_directories(geo PUBLIC ${gwb_global_dir}/include/Geo)
target_link_libraries(geo Eigen3::Eigen)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
# source_group_by_dir_prefixed(src SRC_FILES)
