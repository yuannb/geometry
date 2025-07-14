
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
add_compile_options(-fPIC)


file(GLOB INC_FILES ${ng_global_dir}/include/*.hpp)
file(GLOB SRC_FILES *.cpp)

message(${INC_FILES})
list(APPEND EXTRA_LIBS ${pre_ng_global_dir}/3rd/MySQL/lib/libmysql.lib)
list(APPEND EXTRA_LIBS log)

list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/3rd/CDT/CDT/include)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/3rd/eigen/include/eigen3)
list(APPEND EXTRA_INCLUDES ${ng_global_dir}/include)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/geometry/include/NurbsTools)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/LoggerManager)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/geometry/include/Geo)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/geometry)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/3rd/libhv/include/hv)
list(APPEND EXTRA_INCLUDES ${pre_ng_global_dir}/3rd/MySQL/include)
message("mysqldir: ${pre_ng_global_dir}3rd/MySQL/include")

add_library(node_graph ${INC_FILES} ${SRC_FILES})
target_link_libraries(node_graph PUBLIC ${EXTRA_LIBS})
target_include_directories(node_graph PUBLIC "${EXTRA_INCLUDES}")
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/3rd/CDT/CDT/include)
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/3rd/eigen/include/eigen3)
# target_include_directories(node_graph PUBLIC ${ng_global_dir}/include)
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/geometry/include/NurbsTools)
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/geometry/include/Geo)
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/3rd/libhv/include/hv)
# target_include_directories(node_graph PUBLIC ${pre_ng_global_dir}/3rd/MySQL/include)


if(MSVC)
    target_compile_options(node_graph PRIVATE /utf-8)
endif()

if(WIN32)
    set(MYSQL_RUNTIME_LIB "libmysql.dll")
else()
    set(MYSQL_RUNTIME_LIB "libmysql.so")
endif()

find_file(MYSQL_DLL
    NAMES ${MYSQL_RUNTIME_LIB}
    PATHS ${pre_ng_global_dir}/3rd/MySQL/lib  # MySQL 库目录（需预先设置）
    REQUIRED
)

add_custom_command(
    TARGET node_graph
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${MYSQL_DLL}"
        $<TARGET_FILE_DIR:node_graph>
    COMMENT "Copying MySQL runtime library..."
)

source_group(include FILES ${INC_FILES})
source_group(src FILES ${SRC_FILES})
