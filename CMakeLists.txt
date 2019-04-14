cmake_minimum_required(VERSION 3.2)
project(ram)

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
endif()
# if (NOT TARGET thread_pool)
    # add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
# endif()
if (NOT TARGET logger)
    add_subdirectory(vendor/logger EXCLUDE_FROM_ALL)
endif()

add_executable(ram
    src/aligner.cpp
    src/main.cpp
    src/minimizers.cpp)

target_link_libraries(ram bioparser logger)

install(TARGETS ram DESTINATION ${CMAKE_INSTALL_BINDIR})
