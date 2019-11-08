cmake_minimum_required(VERSION 2.9)
project(Demo)

include_directories("${PROJECT_SOURCE_DIR}/include")

set(CMAKE_BUILD_TYPE Release)

ADD_DEFINITIONS(-DDEBUG)

#MESSAGE(STATUS "the subdir is ${CMAKE_SOURCE_DIR}/lib/")
aux_source_directory("${CMAKE_SOURCE_DIR}/lib" DIR_LIB_SRCS)
#FILE(GLOB DIR_LIB_SRCS  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "./lib/*.c")
#MESSAGE(STATUS "the files in subdir is ${DIR_LIB_SRCS}")
add_library(mylibs SHARED ${DIR_LIB_SRCS})


file(GLOB mains RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/ *.c)

foreach(mainf IN LISTS mains)
    string(REPLACE ".c" "" mainexe ${mainf})
    add_executable(${mainexe} ${mainf})
    target_link_libraries(${mainexe} mylibs -lm) 
endforeach(mainf)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}) 

file(COPY plot.py DESTINATION ${PROJECT_BINARY_DIR})