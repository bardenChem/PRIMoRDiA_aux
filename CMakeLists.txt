cmake_minimum_required(VERSION 3.0)
project(primordia_aux)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(./include)

# Find chemfiles in /usr/local
find_path(CHEMFILES_INCLUDE_DIR chemfiles.hpp
    PATHS /usr/local/include /usr/include
    REQUIRED
)

find_library(CHEMFILES_LIBRARY
    NAMES chemfiles libchemfiles
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64
    REQUIRED
)

message(STATUS "Found chemfiles includes: ${CHEMFILES_INCLUDE_DIR}")
message(STATUS "Found chemfiles library: ${CHEMFILES_LIBRARY}")

# Source files
file(GLOB SOURCES "src/*.cpp")

# Create executable
add_executable(primordia_aux ${SOURCES})

# Include directories
target_include_directories(primordia_aux PRIVATE 
    ${CHEMFILES_INCLUDE_DIR}
)

# Compiler options
target_compile_options(primordia_aux PRIVATE 
    -fopenmp 
    -Wall
    -Wextra
    $<$<CONFIG:Release>:-O3 -DNDEBUG >
    $<$<CONFIG:Debug>:-g -O0 -pg>
)

# Link libraries
target_link_libraries(primordia_aux PRIVATE
    ${CHEMFILES_LIBRARY}
    stdc++fs
)

# Link options
target_link_options(primordia_aux PRIVATE -fopenmp)
