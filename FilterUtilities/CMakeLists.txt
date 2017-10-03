include_directories (${PROJECT_SOURCE_DIR}/FilterUtilities)
set(FILTER_UTILITIES_SOURCES FilterUtilities.cpp GolaySavitzkyCoeff.cpp FFTCutoff.cpp)
add_library(FilterUtilities ${FILTER_UTILITIES_SOURCES})