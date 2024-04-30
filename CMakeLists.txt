project(CGDI-PROJECT)

CMAKE_MINIMUM_REQUIRED(VERSION 3.13)

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
# Geometry processing

include (geometry-central)
include (polyscope)
add_subdirectory(src)

