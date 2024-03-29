cmake_minimum_required(VERSION 3.0.0)
project(dDimensionalSimulations)

add_definitions(-DDIMENSION=9)

set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -DSINGLETHREADED -Xpreprocessor -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Xpreprocessor -fopenmp")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS}
                                        -arch=sm_30
                                        -gencode=arch=compute_30,code=sm_30
                                        -gencode=arch=compute_50,code=sm_50
                                        -gencode=arch=compute_52,code=sm_52
                                        -gencode=arch=compute_60,code=sm_60
                                        -gencode=arch=compute_61,code=sm_61
                                        -gencode=arch=compute_61,code=compute_61
                                        ")
                                        #")#-gencode arch=compute_30,code=sm_30 -Wno-deprecated-gpu-targets ")#-Xptxas -O3 -Xcompiler -O3")

find_package(CUDA REQUIRED)

if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -lineinfo -Xptxas --generate-line-info")
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

message("cxx flags = " "${CMAKE_CXX_FLAGS}")
message("nvcc flags = " "${CUDA_NVCC_FLAGS}")
message("lib dir = " "${CUDA_LIBRARY_DIRS}")
message("culibs = " "${CUDA_LIBRARIES}")
#message("nvcc flags are  " ${CUDA_NVCC_FLAGS})
#message("cuda root is " ${CUDA_INCLUDE_DIRS})

add_subdirectory(src)

include_directories(inc ${CUDA_INCLUDE_DIRS} src src/utilities src/simulation src/model src/forces src/updaters)
link_directories(${CUDA_LIBRARY_DIRS})

set(myLibs
    force
    forceGPU
    updater
    updaterGPU
    simulation
    model
    modelGPU
    utility
    utilityGPU
    )
message("libraries = " "${myLibs}")

add_executable(simulationExecutable.out dDimensionalSimulation.cpp)
if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    target_link_libraries(simulationExecutable.out ${myLibs} ${CUDA_LIBRARIES} "/usr/local/cuda/lib64/libnvToolsExt.so")
else()
target_link_libraries(simulationExecutable.out ${myLibs} ${CUDA_LIBRARIES} gomp)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")
