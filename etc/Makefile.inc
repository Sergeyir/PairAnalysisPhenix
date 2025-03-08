CXX=g++
CXX_COMMON_LIB=-Wall -Werror -Wpedantic -pipe -O2 -fPIC -c
CXX_COMMON_EXE=-Wall -Wpedantic -pipe -O2 -g -ldl

JSON_INCLUDE=-I /usr/include/jsoncpp
JSON_LIB=-ljsoncpp

YAML_CPP_INCLUDE=-I ./yaml-cpp/include
YAML_CPP_LIB=-L ./yaml-cpp/build -Wl,-rpath,./yaml-cpp/build -lyaml-cpp

CPP_TOOLS_LIB_DIR=./CppTools/lib
CPP_TOOLS_SRC_DIR=./CppTools/src
CPP_TOOLS_INCLUDE=-I ./CppTools/include
CPP_TOOLS_LIB=-L $(CPP_TOOLS_LIB_DIR) -Wl,-rpath,$(CPP_TOOLS_LIB_DIR) -lBox -lTime -lTable -lIOTools -lStrTools -lErrorHandler 

ROOT_TOOLS_LIB_DIR=./ROOTTools/lib
ROOT_TOOLS_SRC_DIR=./ROOTTools/src
ROOT_TOOLS_INCLUDE=-I ./ROOTTools/include
ROOT_TOOLS_LIB=-L $(ROOT_TOOLS_LIB_DIR) -Wl,-rpath,$(ROOT_TOOLS_LIB_DIR) -lTCanvasPrinter -lThrObj -lGUIDistrCutter2D

PBAR_LIB_DIR=./ProgressBar/lib
PBAR_SRC_DIR=./ProgressBar/src
PBAR_INCLUDE=-I ./ProgressBar/include
PBAR_LIB=-L $(PBAR_LIB_DIR) -Wl,-rpath,$(PBAR_LIB_DIR) -lPBar

ROOT_INCLUDE=-I${ROOT_PATH}/include
ROOT_CONFIG=${ROOT_PATH}/bin/root-config

PAIR_ANALYSIS_INCLUDE=-I ./include
PAIR_ANALYSIS_LIB=-L ./lib -Wl,-rpath,./lib -lEffTreeReader -lEmbTreeReader -lSTrackFun -lPTrackFun -lIdentFun -lInputReader -lDataMethodsSelector

ALL_INCLUDE=$(JSON_INCLUDE) $(YAML_CPP_INCLUDE) $(CPP_TOOLS_INCLUDE) $(ROOT_INCLUDE) $(PBAR_INCLUDE) $(ROOT_TOOLS_INCLUDE) $(PAIR_ANALYSIS_INCLUDE)
ALL_LIB=$(JSON_LIB) $(YAML_CPP_LIB) $(CPP_TOOLS_LIB) `$(ROOT_CONFIG) --glibs` $(PBAR_LIB) $(ROOT_TOOLS_LIB) $(PAIR_ANALYSIS_LIB)
