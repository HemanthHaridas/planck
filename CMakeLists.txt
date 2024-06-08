cmake_minimum_required(VERSION 3.10)
project(planck_cpp)

# Set the minimum C++ standard and use the most recent standard supported by the compiler
set(CMAKE_CXX_STANDARD 17)  # Minimum C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)  # Disable compiler-specific extensions

# Print the C++ standard being used
message(STATUS "Using C++ Standard: ${CMAKE_CXX_STANDARD}")

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ftree-vectorize -Wall")

# Gather source files
file(GLOB BASE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB LOOKUP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup/*.cpp)
file(GLOB ONE_ELECTRON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/one_electron/*.cpp)
file(GLOB TWO_ELECTRON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/two_electron/*.cpp)
file(GLOB MATH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/math/*.cpp)

# Add the executable target
add_executable(planck 
    src/planck_main.cpp
    ${BASE_SOURCES}
    ${LOOKUP_SOURCES}
    ${ONE_ELECTRON_SOURCES}
    ${TWO_ELECTRON_SOURCES}
    ${MATH_SOURCES}
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
