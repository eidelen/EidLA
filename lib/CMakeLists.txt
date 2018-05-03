cmake_minimum_required(VERSION 3.0)

PROJECT(mylinalglib)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wpedantic -msse4")

INCLUDE_DIRECTORIES(inc)
FILE(GLOB_RECURSE  MYLINALG_INCLUDES        inc/*.hpp)
FILE(GLOB_RECURSE  MYLINALG_INCLUDES        inc/*.h)

option(TESTMYLINALG  "TEST" OFF)
IF(${TESTMYLINALG})
    MESSAGE(STATUS "Tests activated")
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})

    FILE(GLOB_RECURSE  MYLINALG_TESTS_INC       test/*.hpp)
    FILE(GLOB_RECURSE  MYLINALG_TESTS_SRC       test/*.cpp)

    add_executable(runTests ${MYLINALG_TESTS_INC} ${MYLINALG_TESTS_SRC})
    target_link_libraries(runTests ${GTEST_LIBRARIES} )
    target_compile_features( runTests PRIVATE cxx_range_for )
    target_compile_options( runTests PRIVATE "-Wpedantic" )
ENDIF()






