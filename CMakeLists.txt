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

option(ram_build_executable "Build ram standalone tool" OFF)

add_library(ram
    src/minimizer_engine.cpp)

if (NOT TARGET thread_pool)
    add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
endif()

target_link_libraries(ram thread_pool)

target_include_directories(ram PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

install(TARGETS ram DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ram DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if (ram_build_executable)
    add_executable(ram_bin
        src/main.cpp)

    if (NOT TARGET bioparser)
        add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
    endif()
    if (NOT TARGET logger)
        add_subdirectory(vendor/logger EXCLUDE_FROM_ALL)
    endif()

    target_link_libraries(ram_bin ram bioparser thread_pool logger)
    set_target_properties(ram_bin PROPERTIES OUTPUT_NAME ram)

    install(TARGETS ram DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(TARGETS ram_bin DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
