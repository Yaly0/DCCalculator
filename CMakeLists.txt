cmake_minimum_required(VERSION 3.12)
project(DCCalculator)

set(CMAKE_CXX_STANDARD 14)

#find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Eigen3)

add_executable(DCCalculator main.cpp)

target_link_libraries(DCCalculator Eigen3::Eigen)
