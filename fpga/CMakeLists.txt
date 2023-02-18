cmake_minimum_required(VERSION 3.3)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../XilinxCMake")

project(MD++ C CXX)

add_executable(md main.cpp)

find_package(VitisHLS)

MESSAGE(STATUS "${VITIS_HLS_INCLUDE_DIR}")
target_include_directories(md PRIVATE ${VITIS_HLS_INCLUDE_DIR})

set_property(TARGET md PROPERTY CXX_STANDARD 14)
set_property(TARGET md PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET md PROPERTY CXX_EXTENSIONS OFF)
