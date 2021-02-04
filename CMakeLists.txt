cmake_minimum_required(VERSION 3.11)

project(ram VERSION 2.0.0
            LANGUAGES CXX
            DESCRIPTION "Ram is a c++ implementation of [minimap](https://github.com/lh3/minimap) with few modifications.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(ram_main_project ON)
endif ()
option(ram_install "Generate install target" OFF)
option(ram_build_exe "Build executable" ${ram_main_project})
option(ram_build_tests "Build unit tests" ${ram_main_project})

if (ram_install)
  find_package(biosoup 0.10.0 REQUIRED)
  find_package(thread_pool 3.0.3 REQUIRED)
else ()
  find_package(biosoup 0.10.0 QUIET)
  find_package(thread_pool 3.0.3 QUIET)

  if (NOT biosoup_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      biosoup
      GIT_REPOSITORY https://github.com/rvaser/biosoup
      GIT_TAG 0.10.0)

    FetchContent_GetProperties(biosoup)
    if (NOT biosoup_POPULATED)
      FetchContent_Populate(biosoup)
      add_subdirectory(
        ${biosoup_SOURCE_DIR}
        ${biosoup_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()

  if (NOT thread_pool_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      thread_pool
      GIT_REPOSITORY https://github.com/rvaser/thread_pool
      GIT_TAG 3.0.3)

    FetchContent_GetProperties(thread_pool)
    if (NOT thread_pool_POPULATED)
      FetchContent_Populate(thread_pool)
      add_subdirectory(
        ${thread_pool_SOURCE_DIR}
        ${thread_pool_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()
endif ()

if (ram_build_exe OR ram_build_tests)
  find_package(bioparser 3.0.13 QUIET)
  if (NOT bioparser_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      bioparser
      GIT_REPOSITORY https://github.com/rvaser/bioparser
      GIT_TAG 3.0.13)

    FetchContent_GetProperties(bioparser)
    if (NOT bioparser_POPULATED)
      FetchContent_Populate(bioparser)
      add_subdirectory(
        ${bioparser_SOURCE_DIR}
        ${bioparser_BINARY_DIR}
        EXCLUDE_FROM_ALL)
    endif ()
  endif ()
endif ()

if (ram_build_tests)
  find_package(GTest 1.10.0 QUIET)
  if (NOT GTest_FOUND)
    include(FetchContent)

    FetchContent_Declare(
      googletest
      GIT_REPOSITORY https://github.com/google/googletest
      GIT_TAG release-1.10.0)

    FetchContent_GetProperties(googletest)
    if (NOT googletest_POPULATED)
      FetchContent_Populate(googletest)
      add_subdirectory(
        ${googletest_SOURCE_DIR}
        ${googletest_BINARY_DIR}
        EXCLUDE_FROM_ALL)
      add_library(GTest::Main ALIAS gtest_main)
    endif ()
  endif ()
endif ()

add_library(ram
  src/minimizer_engine.cpp)
add_library(ram::ram ALIAS ram)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(ram
  biosoup::biosoup
  thread_pool::thread_pool)

if (ram_install)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(
    TARGETS ram
    EXPORT ${PROJECT_NAME}Targets)
  install(
    DIRECTORY include/ram
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  install(
    EXPORT ${PROJECT_NAME}Targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
  install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
          ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
endif ()

if (ram_build_exe)
  add_executable(ram_exe
    src/main.cpp)

  target_link_libraries(ram_exe
    ram
    bioparser::bioparser)

  target_compile_definitions(ram_exe PRIVATE VERSION="${PROJECT_VERSION}")
  set_property(TARGET ram_exe PROPERTY OUTPUT_NAME ram)

  install(TARGETS ram_exe DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (ram_build_tests)
  add_executable(ram_test
    test/minimizer_engine_test.cpp)

  target_link_libraries(ram_test
    ram
    bioparser::bioparser
    GTest::Main)

  target_compile_definitions(ram_test
    PRIVATE TEST_DATA="${PROJECT_SOURCE_DIR}/test/data/sample.fasta.gz")
endif ()
