cmake_minimum_required(VERSION 3.9)
project(TallBoAnalysis)

set(CMAKE_CXX_STANDARD 11)

# root
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories (${PROJECT_SOURCE_DIR}/FilterUtilities)
include_directories (${PROJECT_SOURCE_DIR}/ROOTUtilities)
add_subdirectory(FilterUtilities)
add_subdirectory(ROOTUtilities)

set(EXTRA_PACKAGES ${EXTRA_PACKAGES} FilterUtilities)
set(EXTRA_PACKAGES ${EXTRA_PACKAGES} ROOTUtilities)
set(SOURCE_FILES main.cpp)
add_executable(TallBoAnalysis ${SOURCE_FILES})

# armadillo by brew (no external link/include needed)
target_link_libraries(TallBoAnalysis armadillo ${ROOT_LIBRARIES} Core ${EXTRA_PACKAGES})