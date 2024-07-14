CXX=g++
CXX_FLAGS_LIB=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c

CPP_TOOLS_INCLUDE=-I ./CppTools/include

ERROR_HANDLER_LIB=-L./lib/ -lErrorHandler
STR_TOOLS_LIB=-L./lib/ -lStrTools

PBAR_LIB=-L./lib -lPBar
PBAR_INCLUDE=-I ./ProgressBar/include
