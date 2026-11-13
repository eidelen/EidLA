cmake_minimum_required(VERSION 3.28)

PROJECT(svd_compression)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(EidLA PATHS ../..)
if(NOT EIDLA_FOUND)
    message(FATAL_ERROR "EidLA could not be found. It is available under https://github.com/eidelen/EidLA")
endif(NOT EIDLA_FOUND)

find_package( OpenCV REQUIRED )
add_definitions(-DOPENCVEIDLA) # to tell eidla to activate opencv export and import

include_directories(${EIDLA_INCLUDE_DIR})
add_executable(svd_compression main.cpp)
target_link_libraries(svd_compression ${OpenCV_LIBS})
target_compile_features(svd_compression PRIVATE cxx_std_17)
