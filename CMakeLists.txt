# meta info
cmake_minimum_required(VERSION 3.9)
project(TallBoAnalysis)
set(CMAKE_CXX_STANDARD 11)

# root package
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# extra packages (written by Dung)
include_directories(${PROJECT_SOURCE_DIR})
add_subdirectory(FilterUtilities)
add_subdirectory(ROOTUtilities)
set(EXTRA_PACKAGES ${EXTRA_PACKAGES} FilterUtilities)
set(EXTRA_PACKAGES ${EXTRA_PACKAGES} ROOTUtilities)

# set main propram
set(SOURCE_FILES main.cpp)
add_executable(TallBoAnalysis ${SOURCE_FILES})

# root (lib + core)
# gsl (standalone)
# fftw
# extra package
target_link_libraries(TallBoAnalysis ${ROOT_LIBRARIES} Core)
target_link_libraries(TallBoAnalysis gsl gslcblas)
target_link_libraries(TallBoAnalysis ${EXTRA_PACKAGES})
target_link_libraries(TallBoAnalysis fftw3 m)