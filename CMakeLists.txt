cmake_minimum_required(VERSION 3.5.0)
project(seqt VERSION 0.1.0 LANGUAGES C CXX)

set(CXX_STANDARD 20)
add_compile_options(--std=c++20)

include_directories(include)

include(CTest)
enable_testing()

add_executable(test_splay test/test_splay.cpp)
add_executable(seqt test/test_seqt.cpp src/seqt.cpp)

add_test(NAME test_splay COMMAND test_splay)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
