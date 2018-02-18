cmake_minimum_required(VERSION 3.9)
project(distance)

set(CMAKE_CXX_STANDARD 14)

include_directories(external)

add_executable(distance main.cpp distance.cpp distance.h util.cpp util.h)
add_executable(run_tests test.cpp test_distance.cpp distance.cpp distance.h util.cpp util.h)

# Needed for Linux; should probably ignore for Windows
target_link_libraries(distance pthread)
target_link_libraries(run_tests pthread)