CXX=g++
CXX_COMMON=-Wpedantic -Wall -O2 -pipe -ldl -g

ROOT_CONFIG=${ROOT}/bin/root-config
ROOT_BIN=${ROOT}/bin
ROOT_INCLUDE=${ROOT}/include
ROOT_LIB=-Wl,-rpath,${ROOT}/lib -L${ROOT}/lib -lCore

PBAR_LIB=-L../../ProgressBar/lib -lPBar
PBAR_INCLUDE=-I ../../ProgressBar/include
