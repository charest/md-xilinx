cmake_minimum_required(VERSION 3.11)
project(MD++ C CXX)

add_executable(md)

set_property(TARGET md PROPERTY CXX_STANDARD 14)
set_property(TARGET md PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET md PROPERTY CXX_EXTENSIONS OFF)

add_subdirectory(host)
