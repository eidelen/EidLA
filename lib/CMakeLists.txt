cmake_minimum_required(VERSION 3.28)

PROJECT(eidlalib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -pthread -Wpedantic -msse4")

INCLUDE_DIRECTORIES(inc inc/pattern)
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
        
    # Download and build gtest
    include(FetchContent)
    FetchContent_Declare(
        googletest
        # Pin the latest release and update it regularly (requires C++17).
        URL https://github.com/google/googletest/archive/refs/tags/v1.18.0.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    FILE(GLOB_RECURSE  MYLINALG_TESTS_INC       test/*.hpp)
    FILE(GLOB_RECURSE  MYLINALG_TESTS_SRC       test/*.cpp)

    add_executable(runTests ${MYLINALG_TESTS_INC} ${MYLINALG_TESTS_SRC}
            test/t_numerical_accuracy.cpp)
    target_link_libraries(runTests ${eidlalibs} gtest_main )
    target_compile_features( runTests PRIVATE cxx_std_17 )
    target_compile_options( runTests PRIVATE "-Wpedantic" )
ENDIF()



