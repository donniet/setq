cmake_minimum_required(VERSION 3.5.0)
project(seq VERSION 0.1.0 LANGUAGES C CXX)

set(CXX_STANDARD 20)
add_compile_options(--std=c++20)

include(CTest)
enable_testing()

add_executable(splay splay.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
