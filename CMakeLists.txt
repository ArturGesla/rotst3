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
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/UMFPACK/Include")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/AMD/Include")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../SuiteSparse-5.12.0/SuiteSparse_config")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/arpackpp-for-eigen/include")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/arpackpp-for-eigen/examples/matrices/nonsym")
list(APPEND INC_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/external/arpackpp-for-eigen/examples/matrices/complex")
 
target_include_directories(rotst PUBLIC ${INC_DIRS})

# # lapack search
find_library(LAPACK_DIR NAMES openblas) #lapack is built in
message (STATUS "lapack found at: ${LAPACK_DIR}")
target_link_libraries(rotst ${LAPACK_DIR})

# superlu search
find_library(SUPERLU_DIR NAMES libsuperlu.so libsuperlu.so.5.3.0 PATHS "/people/gesla/git/superlu/build/SRC")
#find_library(SUPERLU_DIR NAMES libsuperlu.so libsuperlu.so.5 libsuperlu.so.5.3.0 PATHS "/people/gesla/git/superlu/build/SRC")
message (STATUS "superlu found at: ${SUPERLU_DIR}")
target_link_libraries(rotst ${SUPERLU_DIR})

# arpack search (probably in some local dir)
find_library(ARPACK_DIR NAMES arpack libarpack.so.2 PATHS ${LOC_LIB_DIRS})
message (STATUS "arpack found at: ${ARPACK_DIR}")
target_link_libraries(rotst ${ARPACK_DIR})

#umfpack search
find_library(UMFPACK_DIR NAMES libumfpack.so libumfpack.so.5 libumfpack.so.5.7.1 PATHS ${LOC_LIB_DIRS})
message (STATUS "umfpack found at: ${UMFPACK_DIR}")
target_link_libraries(rotst ${UMFPACK_DIR})
