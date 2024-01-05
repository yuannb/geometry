
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_compile_options(-fPIC)

file(GLOB INC_FILES ${gwb_global_dir}/include/Geo/*.h)
file(GLOB SRC_FILES *.cpp)

add_library(geo ${INC_FILES} ${SRC_FILES})

target_include_directories(geo PUBLIC ${gwb_global_dir}/include/Geo)
target_include_directories(geo PUBLIC ${gwb_global_dir}/include/NurbsTools)
target_include_directories(geo PUBLIC ${pre_gwb_global_dir}/3rd/eigen/include/eigen3)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
# source_group_by_dir_prefixed(src SRC_FILES)
