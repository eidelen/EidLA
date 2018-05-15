cmake_minimum_required(VERSION 3.0)

PROJECT(eidlalib)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wpedantic -msse4")

INCLUDE_DIRECTORIES(inc)
FILE(GLOB_RECURSE  MYLINALG_INCLUDES        inc/*.hpp)
FILE(GLOB_RECURSE  MYLINALG_INCLUDES        inc/*.h)

set(eidlalibs )

option(USEOPENCV  "OpenCV compatible" OFF)
IF(${USEOPENCV})
    find_package( OpenCV REQUIRED )
    add_definitions(-DOPENCVEIDLA)
    set(eidlalibs ${eidlalibs} ${OpenCV_LIBS} )
    MESSAGE(STATUS "Executables required to link OpenCV!")
ENDIF(${USEOPENCV})

option(TESTMYLINALG  "TEST" OFF)
IF(${TESTMYLINALG})
    MESSAGE(STATUS "Tests activated")
    find_package(GTest REQUIRED)
    include_directories(${GTEST_INCLUDE_DIRS})

    FILE(GLOB_RECURSE  MYLINALG_TESTS_INC       test/*.hpp)
    FILE(GLOB_RECURSE  MYLINALG_TESTS_SRC       test/*.cpp)

    add_executable(runTests ${MYLINALG_TESTS_INC} ${MYLINALG_TESTS_SRC})
    target_link_libraries(runTests ${GTEST_LIBRARIES} ${eidlalibs} )
    target_compile_features( runTests PRIVATE cxx_range_for )
    target_compile_options( runTests PRIVATE "-Wpedantic" )
ENDIF()






