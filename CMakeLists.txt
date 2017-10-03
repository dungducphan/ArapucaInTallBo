cmake_minimum_required(VERSION 3.9)
project(TallBoAnalysis)

set(CMAKE_CXX_STANDARD 11)

# root
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

set(SOURCE_FILES main.cpp ROOTUtilities/Analysis.h ROOTUtilities/Analysis.cpp FilterUtilities/GolaySavitzkyCoeff.cpp FilterUtilities/GolaySavitzkyCoeff.h FilterUtilities/FilterUtilities.cpp FilterUtilities/FilterUtilities.h)
add_executable(TallBoAnalysis ${SOURCE_FILES})

# armadillo by brew (no external link/include needed)
target_link_libraries(TallBoAnalysis armadillo ${ROOT_LIBRARIES} Core)