# if(MSVC)
#     set(CMAKE_GENERATOR_PLATFORM "x64" CACHE STRING "Target Platform" FORCE)
# endif()
cmake_minimum_required(VERSION 3.10)
project(colmap2pba LANGUAGES CXX)

# C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# ----------------------------
# OpenCV
# ----------------------------
find_package(OpenCV REQUIRED)
# ----------------------------
# Eigen
# ----------------------------
set(3RDPARTY_DIR ${CMAKE_SOURCE_DIR}/3rdparty)
message("3RDPARTY_DIR = ${3RDPARTY_DIR}")
# Eigen 设置
set(EIGEN3_INCLUDE_DIR "${3RDPARTY_DIR}/eigen3")


# ----------------------------
# 源文件
# ----------------------------
set(SOURCES
    "src/mycommon.cpp"
    "src/myconvert.cpp"
    "src/myio.cpp"
    "src/triangulation.cpp"
    "src/colmap2pba.cpp"
)

# 创建可执行程序
add_executable(colmap2pba ${SOURCES})

set(COLMAP2PBA_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# 头文件路径
target_include_directories(colmap2pba PRIVATE
    # ${CMAKE_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
    # ${OpenCV_INCLUDE_DIRS}
    ${COLMAP2PBA_INCLUDE_DIR}
    # ${3RDPARTY_DIR}
)
target_link_libraries(colmap2pba PRIVATE ${OpenCV_LIBS})

# add_custom_command(TARGET colmap2pba POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#     "${OpenCV_DIR}/x64/vc16/bin/opencv_world4110.dll"
#     $<TARGET_FILE_DIR:colmap2pba>)

# if (WIN32 AND OpenCV_RUNTIME_DLLS)
#     add_custom_command(TARGET colmap2pba POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         ${OpenCV_RUNTIME_DLLS}
#         $<TARGET_FILE_DIR:colmap2pba>
#     )
# endif()