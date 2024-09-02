CXX=g++
CXX_FLAGS=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c
CXX_COMMON=-Wall -Wpedantic -pipe -O2 -ldl -g

CPP_TOOLS_INCLUDE=-I ./CppTools/include
CPP_TOOLS_LIB=-L ./lib -lErrorHandler -lStrTools -lIOTools -lBox -lTime -lTable

ROOT_TOOLS_INCLUDE=-I ./ROOTTools/include
ROOT_TOOLS_LIB=-L ./ROOTTools/lib -l

PBAR_INCLUDE=-I ./ProgressBar/include
PBAR_LIB=-L ./ProgressBar/lib -lPBar

SIM_ANALYSIS_INCLUDE=-I ./SimAnalysis/include
SIM_ANALYSIS_LIB=-L ./SimAnalysis/lib

ROOT_INCLUDE=${ROOT}/include
ROOT_LIB=-Wl,-rpath,${ROOT}/lib -L${ROOT}/lib -lCore
ROOT_CONFIG=${ROOT}/bin/root-config
