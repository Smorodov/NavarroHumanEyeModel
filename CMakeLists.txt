cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

project(Human-Eye-Model-Ray-Tracing)

    FIND_PACKAGE(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()

include_directories(gdt)
include_directories(gdt/gdt)
include_directories(gdt/gdt/math)
include_directories(gdt/gdt/random)
add_subdirectory(gdt)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

set(SOURCES
src/eye.cpp
src/main.cpp
)

set(HEADERS
src/eye.h
)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)
SOURCE_GROUP("Source Files" FILES ${SOURCES})
SOURCE_GROUP("Header Files" FILES ${HEADERS})
# -----------------------------
#
# -----------------------------
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} gdt)
# -----------------------------
set_target_properties( ${PROJECT_NAME} PROPERTIES  VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROJECT_NAME}>")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_SOURCE_DIR}/build/Debug)
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_SOURCE_DIR}/build/Release)
