cmake_minimum_required(VERSION 3.10)

# MODIFY: Project name
project(criticalPoint)

# MODIFY: This needs to point to the Lib directory
set(PATH_TO_DFT_LIB ${PROJECT_SOURCE_DIR}/../../Lib)

# MODIFY: this is the app file
add_executable(${PROJECT_NAME}  ${PROJECT_SOURCE_DIR}/criticalPoint.cpp)

#MODIFY: enable these lines if you want to use the WHDF potential
set(POTENTIAL_WHDF true)

#MODIFY: enable these lines if you want to use an analytic evaluation of weights
set(ANALYTIC_WEIGHTS false)

################### End of customizations

## Set up the compilation
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
set(CMAKE_CXX_EXTENSIONS ON)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

## Find other packages needed
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/../../scripts)

find_package(GSL REQUIRED)
find_package(OpenMP)
find_package(Armadillo  REQUIRED)

if(OpenMP_CXX_FOUND)
  target_compile_definitions(${PROJECT_NAME} PUBLIC USE_OMP)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

## Get the classicalDFT library and headers
find_library(CLASSICALDFT classicalDFT  PATHS ${PATH_TO_DFT_LIB})
include_directories(${PATH_TO_DFT_LIB}/include)

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CLASSICALDFT}
  "gsl"
  "gslcblas"
  "armadillo"
  "fftw3_omp"
  "fftw3"
  )

if(POTENTIAL_WHDF)
  target_compile_definitions(${PROJECT_NAME} PUBLIC POTENTIAL_WHDF)
endif()

if(ANALYTIC_WEIGHTS)
  target_compile_definitions(${PROJECT_NAME} PUBLIC ANALYTIC_WEIGHTS)
endif()


