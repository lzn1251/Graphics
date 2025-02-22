cmake_minimum_required(VERSION 3.21)

project(GRAPHICS)

include_directories(${PROJECT_SOURCE_DIR}/dependency)

add_subdirectory(source/BSDF)

add_executable(Graphics source/main.cpp)

target_link_libraries(Graphics BSDF)