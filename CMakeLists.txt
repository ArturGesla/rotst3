cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#SET(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp")
#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

#set(ENV{ROTST_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/build/")

project(rotst3 VERSION 1.0)

file(GLOB_RECURSE SRC_FILES library/src/*.cpp)

add_executable(rotst drivers/main.cpp ${SRC_FILES})

set(INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/library/inc")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../eigen")
 
target_include_directories(rotst PUBLIC ${INC_DIRS})
