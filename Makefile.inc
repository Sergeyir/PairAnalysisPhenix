CXX=g++
CXX_FLAGS=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c


CPP_TOOLS_INCLUDE=-I ./CppTools/include
ERROR_HANDLER_LIB=-L./lib/ -lErrorHandler
STR_TOOLS_LIB=-L./lib/ -lStrTools
IOTOOLS_LIB=-L./lib/ -lIOTools

ROOT_TOOLS_INCLUDE=-I ./ROOTTools/include
THR_OBJ_LIB=-I ./ROOTTools/include

PBAR_INCLUDE=-I ./ProgressBar/include
PBAR_LIB=-L./lib -lPBar

SIM_ANALYSIS_INCLUDE=-I ./SimAnalysis/include
EMB_TREE_READER_LIB=-L./lib -lEmbTreeReader
S_TRACK_FUN_LIB=-L./lib -lSTrackFun

ROOT_INCLUDE=${ROOT}/include
ROOT_LIB=-Wl,-rpath,${ROOT}/lib -L${ROOT}/lib -lCore
ROOT_CONFIG=${ROOT}/bin/root-config
