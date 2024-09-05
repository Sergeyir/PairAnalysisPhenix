# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

all: all_libs AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack

all_libs: ErrorHandler StrTools Time IOTools Box Table ProgressBar TCanvasPrinter GUIFit EffTreeReader EmbTreeReader STrackFun

# CppTools

CppToolsLib: 
	mkdir -p CppTools/lib

ErrorHandler: CppTools/src/ErrorHandler.cpp CppToolsLib
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o
	
StrTools: CppTools/src/StrTools.cpp CppToolsLib
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

IOTools: CppTools/src/IOTools.cpp ErrorHandler StrTools
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Box: CppTools/src/Box.cpp IOTools
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools -lIOTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Table: CppTools/src/Table.cpp IOTools
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools -lIOTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Time: CppTools/src/Time.cpp CppToolsLib
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

# ProgressBar

ProgressBarLib:
	mkdir -p ProgressBar/lib

ProgressBar: ProgressBar/src/PBar.cpp ProgressBarLib
	$(CXX) ProgressBar/src/PBar.cpp $(CXX_COMMON_LIB) \
	-o ProgressBar/lib/PBar.o && \
	ar rcs ProgressBar/lib/libPBar.a ProgressBar/lib/PBar.o && \
	$(CXX) -shared -o ProgressBar/lib/PBar.so ProgressBar/lib/PBar.o

# ROOTTools

ROOTToolsLib:
	mkdir -p ROOTTools/lib

TCanvasPrinter: ROOTTools/src/TCanvasPrinter.cpp ROOTToolsLib
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

GUIFit: ROOTTools/src/GUIFit.cpp ErrorHandler IOTools ROOTToolsLib
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) -L./CppTools/lib -l ErrorHandler -lIOTools \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

# SimAnalysis

SimAnalysisLib:
	mkdir -p SimAnalysis/lib

EffTreeReader: SimAnalysis/src/EffTreeReader.cpp SimAnalysisLib
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o

EmbTreeReader: SimAnalysis/src/EmbTreeReader.cpp SimAnalysisLib
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o

STrackFun: SimAnalysis/src/STrackFun.cpp SimAnalysisLib
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o
	
# SimAnalysis executables

AnalyzeEmbedding: SimAnalysis/src/AnalyzeEmbedding.cpp all_libs bin
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

AnalyzeHeatMaps: SimAnalysis/src/AnalyzeHeatMaps.cpp all_libs bin
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

AnalyzeSingleTrack: SimAnalysis/src/AnalyzeSingleTrack.cpp all_libs bin
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

# other

bin:
	mkdir bin

clean: 
	rm -r bin ; \
	rm -r CppTools/lib ; \
	rm -r ROOTTools/lib ; \
	rm -r ProgressBar/lib ; \
	rm -r SimAnalysis/lib
