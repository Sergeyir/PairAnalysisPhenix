# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

HeatMapper: src/HeatMapper.cpp
	$(CXX) $@.cpp -o $@.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

EmbAnalyze: src/EmbAnalyze.cpp
	$(CXX) $@.cpp -o $@.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

SingleAnalyze: src/SingleAnalyze.cpp
	$(CXX) $@.cpp -o $@.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

WidthAnalyze: src/WidthAnalyze.cpp
	$(CXX) $@.cpp -o $@.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

PairAnalyze: src/PairAnalyze.cpp
	$(CXX) $@.cpp -o $@.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

all : HeatMapper.cpp EmbAnalyze.cpp SingleAnalyze.cpp WidthAnalyze.cpp PairAnalyze.cpp
	$(CXX) HeatMapper.cpp -o HeatMapper.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` & \
	$(CXX) EmbAnalyze.cpp -o EmbAnalyze.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` & \
	$(CXX) SingleAnalyze.cpp -o SingleAnalyze.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` & \
	$(CXX) WidthAnalyze.cpp -o WidthAnalyze.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` & \
	$(CXX) PairAnalyze.cpp -o PairAnalyze.exe $(CXX_COMMON) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

# Clean.
clean:
	rm *.exe
