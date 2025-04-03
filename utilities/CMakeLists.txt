cmake_minimum_required(VERSION 3.10)
project(planck_utilities)

# Set the minimum C++ standard and use the most recent standard supported by the compiler
set(CMAKE_CXX_STANDARD 17)  # Minimum C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Print the C++ standard being used
message(STATUS "Using C++ Standard: ${CMAKE_CXX_STANDARD}")

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -ftree-vectorize -Wall -pedantic")

# find_package(Boost REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

# Gather source files
file(GLOB MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

# Add the executable target
add_executable(planck_utilities 
    ${MAIN_SOURCES}
)

# Include directories
target_include_directories(planck PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities
)

# Install the executable
install(TARGETS planck_utilities DESTINATION)