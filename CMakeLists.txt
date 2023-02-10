cmake_minimum_required(VERSION 3.3)
project(MD++ C CXX)

add_executable(md main.cpp)

set_property(TARGET md PROPERTY CXX_STANDARD 14)
set_property(TARGET md PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET md PROPERTY CXX_EXTENSIONS OFF)
