cmake_minimum_required (VERSION 3.5)

project (ram-mapper)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++1y")
set (source_dir "${PROJECT_SOURCE_DIR}/src/")
set (CMAKE_BUILD_TYPE Release)

#add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
#add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

#add_library(pink_alignment src/pink_alignment.cpp)
#add_library(pink_minimizers src/pink_minimizers.cpp)

add_executable (ram-mapper src/ram_mapper.cpp)

#target_link_libraries(pink_mapper bioparser pink_alignment pink_minimizers thread_pool pthread)
