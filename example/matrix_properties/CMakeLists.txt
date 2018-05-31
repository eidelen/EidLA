cmake_minimum_required(VERSION 3.0)

PROJECT(eidla_example)

find_package(EidLA PATHS ../..)
if(NOT EIDLA_FOUND)
    message(FATAL_ERROR "EidLA could not be found. It is available under https://github.com/eidelen/EidLA")
endif(NOT EIDLA_FOUND)

include_directories(${EIDLA_INCLUDE_DIR})
add_executable(eidla_example main.cpp)
target_compile_features(eidla_example PRIVATE cxx_range_for)
