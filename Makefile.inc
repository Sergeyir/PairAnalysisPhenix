CXX=g++
CXX_COMMON_LIB=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c
CXX_COMMON_EXE=-Wall -Wpedantic -pipe -O2 -g -ldl

CPP_TOOLS_INCLUDE=-I ./CppTools/include
CPP_TOOLS_LIB=-L./CppTools/lib -lBox -lTime -lTable -lIOTools -lStrTools -lErrorHandler 

ROOT_TOOLS_INCLUDE=-I ./ROOTTools/include
ROOT_TOOLS_LIB=-L ./ROOTTools/lib -lTCanvasPrinter -lGUIFit

PBAR_INCLUDE=-I ./ProgressBar/include
PBAR_LIB=-L ./ProgressBar/lib -lPBar

SIM_ANALYSIS_INCLUDE=-I ./SimAnalysis/include
SIM_ANALYSIS_LIB=-L ./SimAnalysis/lib -lEffTreeReader -lEmbTreeReader -lSTrackFun

ROOT_LIB=-Wl,-rpath,${ROOT}/lib -L${ROOT}/lib -lCore
ROOT_CONFIG=${ROOT}/bin/root-config
