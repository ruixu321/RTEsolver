cmake_minimum_required(VERSION 2.8.9)

#Bring the headers, such as Student.h into the project
include_directories(include)

set (CMAKE_CXX_STANDARD 14)

#Source files
file(GLOB SOURCES "src/*.cpp")

#name of the project
project (radiation)

#add the executables
add_executable(radiation ${SOURCES})
