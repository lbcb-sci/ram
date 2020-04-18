cmake_minimum_required(VERSION 3.9)

project(ram VERSION 1.0.2
            LANGUAGES CXX
            DESCRIPTION "Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT TARGET biosoup)
  add_subdirectory(vendor/bioparser/vendor/biosoup EXCLUDE_FROM_ALL)
endif ()
if (NOT TARGET thread_pool)
  add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
endif ()
add_library(${PROJECT_NAME} src/minimizer_engine.cpp)
target_link_libraries(${PROJECT_NAME} biosoup thread_pool)

target_include_directories(${PROJECT_NAME}
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

option(ram_build_executable "Build ram executable" OFF)
if (ram_build_executable)
  if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
  endif ()
  add_executable(${PROJECT_NAME}_exe src/main.cpp)
  target_link_libraries(${PROJECT_NAME}_exe ${PROJECT_NAME} bioparser)
  target_compile_definitions(${PROJECT_NAME}_exe
    PRIVATE RAM_VERSION="v${PROJECT_VERSION}")
  set_target_properties(${PROJECT_NAME}_exe
    PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

  include(GNUInstallDirs)
  install(TARGETS ${PROJECT_NAME}_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

option(ram_build_tests "Build ram unit tests" OFF)
if (ram_build_tests)
  find_package(GTest REQUIRED)
  if (NOT TARGET bioparser)
    add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
  endif ()
  add_executable(${PROJECT_NAME}_test test/minimizer_engine_test.cpp)
  target_link_libraries(${PROJECT_NAME}_test ${PROJECT_NAME} bioparser GTest::Main)
  target_compile_definitions(${PROJECT_NAME}_test
    PRIVATE RAM_DATA_PATH="${PROJECT_SOURCE_DIR}/test/data/sample.fasta.gz")
endif ()
