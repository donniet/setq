cmake_minimum_required(VERSION 3.5.0)
project(seq VERSION 0.1.0 LANGUAGES C CXX)

set(CXX_STANDARD 20)

include(CTest)
enable_testing()

add_executable(seq main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
