cmake_minimum_required(VERSION 3.0)
project(libnm)

add_subdirectory(nm)

if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    add_subdirectory(tests)
endif()
