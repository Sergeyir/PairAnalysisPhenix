cmake_minimum_required(VERSION 3.0...3.31.6 FATAL_ERROR)

project(PairAnalysisPhenix)

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)


set(CMAKE_CXX_FLAGS "-Wall -Werror -Wpedantic -pipe -O2")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_BUILD_TYPE Debug)

find_package(ROOT QUIET)

if(${ROOT_FOUND})
   message(STATUS "ROOT package found")
   link_libraries(${ROOT_LIBRARIES})
elseif(EXISTS "${ROOT_PATH}/bin/root")
   message(STATUS "ROOT path found")
   set(ROOT_INCLUDE_DIRS "${ROOT_PATH}/include")
   set(ROOT_root_CMD "${ROOT_PATH}/bin/root")
   execute_process(COMMAND ${ROOT_root_CMD}-config --glibs OUTPUT_VARIABLE ROOT_LIB_FLAGS)
   string(STRIP ${ROOT_LIB_FLAGS} ROOT_LIB_FLAGS)
   set(CMAKE_SHARED_LINKER_FLAGS ${ROOT_LIB_FLAGS})
   set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LIB_FLAGS})
else()
   message(FATAL_ERROR "ROOT not found: install it via package manager or set ROOT_PATH if you installed it manually")
endif()

#if (EXISTS /usr/include/json)
#   message(STATUS "jsoncpp package found")
#elseif (EXISTS /usr/include/jsoncpp)
#   message(STATUS "libjsoncpp package found")
#   include_directories(${JSON_CPP_PATH}/include/jsoncpp)
#elseif (EXISTS ${JSON_CPP_PATH}/include/json)
#   message(STATUS "jsoncpp path found")
#   link_directories(${JSON_CPP_PATH}/build/lib)
#   include_directories(${JSON_CPP_PATH}/include)
#else()
#   message(FATAL_ERROR "jsoncpp/libjsoncpp not found: install it via package manager or set environment variable JSON_CPP_PATH if you installed it manually")
#endif()

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/CppTools/include ./ROOTTools/include ${CMAKE_SOURCE_DIR}/ProgressBar/include)

add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)
add_subdirectory(CppTools EXCLUDE_FROM_ALL)
add_subdirectory(ROOTTools)
add_subdirectory(ProgressBar EXCLUDE_FROM_ALL)

#link_libraries(jsoncpp)
link_libraries(ErrorHandler StrTools Time IOTools Box Table)
link_libraries(TCanvasTools ThrObj)
link_libraries(PBar)
link_libraries(yaml-cpp)

add_library(InputYAMLReader ${CMAKE_SOURCE_DIR}/src/InputYAMLReader.cpp)
add_library(SimTreeReader ${CMAKE_SOURCE_DIR}/src/SimTreeReader.cpp)
add_library(SingleTrackFunc ${CMAKE_SOURCE_DIR}/src/SingleTrackFunc.cpp)
add_library(PairTrackFunc ${CMAKE_SOURCE_DIR}/src/PairTrackFunc.cpp)
add_library(DeadMapCutter ${CMAKE_SOURCE_DIR}/src/DeadMapCutter.cpp)
add_library(SimSigmalizedResiduals ${CMAKE_SOURCE_DIR}/src/SimSigmalizedResiduals.cpp)
add_library(SimM2Identificator ${CMAKE_SOURCE_DIR}/src/SimM2Identificator.cpp)
add_library(FitFunc ${CMAKE_SOURCE_DIR}/src/FitFunc.cpp)

link_libraries(InputYAMLReader)

add_executable(AnalyzeSingleTrack ${CMAKE_SOURCE_DIR}/src/AnalyzeSingleTrack.cpp)
add_executable(AnalyzeWidthlessResonance ${CMAKE_SOURCE_DIR}/src/AnalyzeWidthlessResonance.cpp)
add_executable(AnalyzeResonance ${CMAKE_SOURCE_DIR}/src/AnalyzeResonance.cpp)
add_executable(CalibrateSimSigmalizedResiduals ${CMAKE_SOURCE_DIR}/src/CalibrateSimSigmalizedResiduals.cpp)

add_executable(DeadMapSys ${CMAKE_SOURCE_DIR}/src/DeadMapSys.cpp)
add_executable(M2IdentFit ${CMAKE_SOURCE_DIR}/src/M2IdentFit.cpp)

add_executable(EstimateGaussianBroadening ${CMAKE_SOURCE_DIR}/src/EstimateGaussianBroadening.cpp)
add_executable(EstimateSingleTrackEff ${CMAKE_SOURCE_DIR}/src/EstimateSingleTrackEff.cpp)

target_link_libraries(SingleTrackFunc SimTreeReader)
target_link_libraries(PairTrackFunc SimTreeReader)
target_link_libraries(AnalyzeSingleTrack SimTreeReader SingleTrackFunc DeadMapCutter SimSigmalizedResiduals SimM2Identificator)
target_link_libraries(AnalyzeWidthlessResonance SimTreeReader SingleTrackFunc PairTrackFunc DeadMapCutter SimSigmalizedResiduals)
target_link_libraries(AnalyzeResonance SimTreeReader SingleTrackFunc PairTrackFunc DeadMapCutter SimSigmalizedResiduals SimM2Identificator)
target_link_libraries(DeadMapSys DeadMapCutter)
