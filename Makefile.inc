CXX=g++
CXX_FLAGS=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c

CPP_TOOLS_INCLUDE=-I ./CppTools/include

ERROR_HANDLER_LIB=-L./lib/ -lErrorHandler
STR_TOOLS_LIB=-L./lib/ -lStrTools
IOTOOLS_LIB=-L./lib/ -lIOTools

PBAR_LIB=-L./lib -lPBar
PBAR_INCLUDE=-I ./ProgressBar/include

ROOT_CONFIG=${ROOT}/bin/root-config
ROOT_BIN=${ROOT}/bin
ROOT_INCLUDE=${ROOT}/include
ROOT_LIB=-Wl,-rpath,${ROOT}/lib -L${ROOT}/lib -lCore
