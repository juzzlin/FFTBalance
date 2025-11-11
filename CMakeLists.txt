cmake_minimum_required(VERSION 3.20)
project(FFTBalance LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Use pkg-config to find libsndfile
find_package(PkgConfig REQUIRED)
pkg_check_modules(SNDFILE REQUIRED sndfile)

# Use pkg-config to find fftw3
pkg_check_modules(FFTW3 REQUIRED fftw3)

# Include directories
include_directories(
    ${SNDFILE_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
)

# Link directories (pkg-config)
link_directories(
    ${SNDFILE_LIBRARY_DIRS}
    ${FFTW3_LIBRARY_DIRS}
)

set(BINARY_NAME fft_balance) 
add_executable(${BINARY_NAME} main.cpp application.cpp fft_balance.cpp)
target_link_libraries(${BINARY_NAME}
    ${SNDFILE_LIBRARIES}
    ${FFTW3_LIBRARIES}
)
target_compile_options(${BINARY_NAME} PRIVATE -Wall -Wextra -O2)

