cmake_minimum_required(VERSION 3.5.0)
project(seqt VERSION 0.1.0 LANGUAGES C CXX)

set(CXX_STANDARD 20)
add_compile_options(--std=c++20 -g)

include_directories(include)

include(CTest)
enable_testing()

add_executable(test_splay test/test_splay.cpp)
add_executable(test_seqt test/test_seqt.cpp src/seqt.cpp)
add_executable(test_manual test/test_manual.cpp)
add_executable(symbol_workspace src/symbol_workspace.cpp)

add_library(partition src/partition.cpp)

add_executable(test_partition test/test_partition.cpp)
target_link_libraries(test_partition partition)

add_test(NAME test_splay COMMAND test_splay)
add_test(NAME test_seqt COMMAND test_seqt)
add_test(NAME test_manual COMMAND test_manual)
add_test(NAME test_partition COMMAND test_partition)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
