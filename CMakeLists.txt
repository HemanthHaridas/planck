cmake_minimum_required(VERSION 3.10)
project(planck_cpp)

# Set the minimum C++ standard and use the most recent standard supported by the compiler
set(CMAKE_CXX_STANDARD 17)  # Minimum C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/opt/homebrew/Cellar/boost/1.87.0_1/include")
# set(CMAKE_CXX_EXTENSIONS OFF)  # Disable compiler-specific extensions

# Print the C++ standard being used
message(STATUS "Using C++ Standard: ${CMAKE_CXX_STANDARD}")

# Set optimization flags based on OS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Wall -pedantic -ftree-vectorize")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin") # macOS
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Wall -pedantic -Rpass=loop-vectorize")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++${CMAKE_CXX_STANDARD} /O2 /W4 /permissive-")
else()
    message(WARNING "Unknown system; defaulting to generic flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Wall -pedantic")
endif()

message(STATUS "Building on : ${CMAKE_SYSTEM_NAME} using following compiler flags : ${CMAKE_CXX_FLAGS}")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Rpass=loop-vectorize -Wall -pedantic")

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Gather source files
file(GLOB MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB BASE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB LOOKUP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup/*.cpp)
file(GLOB HELPER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/helper/*.cpp)
file(GLOB ONE_ELECTRON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/one_electron/*.cpp)
file(GLOB TWO_ELECTRON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/two_electron/*.cpp)
file(GLOB MATH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/math/*.cpp)

# Add the executable target
add_executable(planck 
    ${MAIN_SOURCES}
    ${BASE_SOURCES}
    ${MATH_SOURCES}
    ${HELPER_SOURCES}
    ${LOOKUP_SOURCES}
    ${ONE_ELECTRON_SOURCES}
    ${TWO_ELECTRON_SOURCES}
)

# Get the dependencies
set(MSYM_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/install)
set(MSYM_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/build)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install)

include(ExternalProject)

ExternalProject_Add(
    libmsym
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym
    GIT_REPOSITORY https://github.com/mcodev31/libmsym.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MSYM_INSTALL_DIR} -DBUILD_SHARED_LIBS:BOOL=OFF -DMSYM_BUILD_EXAMPLES:BOOL=OFF
)

ExternalProject_Get_Property(libmsym install_dir)

# Include directories
target_include_directories(planck PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/helper
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/one_electron
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${MSYM_INSTALL_DIR}/include  # Include the external library headers
)

# Link the external library
target_link_libraries(planck ${MSYM_INSTALL_DIR}/lib/libmsym.a)

# Add dependencies to ensure the external project is built first
add_dependencies(planck libmsym)

# Install the executable
install(TARGETS planck DESTINATION)

# copy the basis set files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/basis DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# write the path of basis sets
message(STATUS "writing planck.defaults file")
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/planck.defaults "${CMAKE_INSTALL_PREFIX}/bin/basis\n")
