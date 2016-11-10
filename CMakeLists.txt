cmake_minimum_required(VERSION 3.2)
project(quantumwalk_cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES ./versions/ultra_quantum_walk/main.cpp)
add_executable(quantumwalk_cpp ${SOURCE_FILES} ./versions/ultra_quantum_walk/UltraQuantumWalk.cpp ./versions/ultra_quantum_walk/UltraQuantumWalk.h)
#set(SOURCE_FILES ./versions/random_medium_semi_inifinite_quantum_walk/main.cpp)
#add_executable(quantumwalk_cpp ${SOURCE_FILES} ./versions/random_medium_semi_inifinite_quantum_walk/QuantumWalk.cpp ./versions/random_medium_semi_inifinite_quantum_walk/QuantumWalk.h)