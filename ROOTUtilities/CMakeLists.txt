include_directories (${PROJECT_SOURCE_DIR}/ROOTUtilities)
set(ROOT_UTILITIES_SOURCES DataAccess.cpp)
add_library(ROOTUtilities ${ROOT_UTILITIES_SOURCES})