cmake_minimum_required(VERSION 3.0)
project("Planck")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# Find Boost components
find_package(Boost REQUIRED COMPONENTS)

# Add source files
file(GLOB EXECUTABLE_SOURCES "src/*.cpp")

# Create an executable for each source file
foreach(EXECUTABLE_SOURCE ${EXECUTABLE_SOURCES})
    get_filename_component(EXECUTABLE_NAME ${EXECUTABLE_SOURCE} NAME_WE)
    add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE})
    # Link libraries
    target_link_libraries(${EXECUTABLE_NAME} Eigen3::Eigen ${Boost_LIBRARIES})
endforeach()

# Set output directory for binary
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)