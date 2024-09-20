# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
I := $(patsubst %/,%,$(dir $(_mkfile_path)))

ifneq ($(words $(MAKECMDGOALS)),1)
.DEFAULT_GOAL = all
%:
	@$(MAKE) $@ --no-print-directory -rRf $(firstword $(MAKEFILE_LIST))
else
ifndef ECHO
T := $(shell $(MAKE) $(MAKECMDGOALS) --no-print-directory \
     -nrRf $(firstword $(MAKEFILE_LIST)) \
     ECHO="COUNTTHIS" | grep -c "COUNTTHIS")
N := x
C = $(words $N)$(eval N := x $N)
ECHO = python3 $(I)/makefile_progress.py --stepno=$C --nsteps=$T
endif

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

.SILENT: CppToolsLib ROOTToolsLib ProgressBarLib SimAnalysisLib ErrorHandler StrTools Time IOTools Box Table ProgressBar TCanvasPrinter GUIFit EffTreeReader EmbTreeReader STrackFun AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack


all: all_libs AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack
	@echo "All done"

all_libs: ErrorHandler StrTools Time IOTools Box Table ProgressBar TCanvasPrinter GUIFit EffTreeReader EmbTreeReader STrackFun

# CppTools

CppToolsLib: 
	@$(ECHO) Creating directory for CppTools libraries
	mkdir -p CppTools/lib

ErrorHandler: CppTools/src/ErrorHandler.cpp CppToolsLib
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o
	
StrTools: CppTools/src/StrTools.cpp CppToolsLib
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

IOTools: CppTools/src/IOTools.cpp ErrorHandler StrTools
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Box: CppTools/src/Box.cpp IOTools
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools -lIOTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Table: CppTools/src/Table.cpp IOTools
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(CPP_TOOLS_INCLUDE) \
	-L ./CppTools/lib -lErrorHandler -lStrTools -lIOTools \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

Time: CppTools/src/Time.cpp CppToolsLib
	@$(ECHO) Compiling $@.cpp from CppTools module
	$(CXX) CppTools/src/$@.cpp $(CXX_COMMON_LIB) \
	-o CppTools/lib/$@.o && \
	ar rcs CppTools/lib/lib$@.a CppTools/lib/$@.o && \
	$(CXX) -shared -o CppTools/lib/$@.so CppTools/lib/$@.o

# ProgressBar

ProgressBarLib:
	@$(ECHO) Creating directory for ProgressBar libraries
	mkdir -p ProgressBar/lib

ProgressBar: ProgressBar/src/PBar.cpp ProgressBarLib
	@$(ECHO) Compiling $@.cpp from ProgressBar module
	$(CXX) ProgressBar/src/PBar.cpp $(CXX_COMMON_LIB) \
	-o ProgressBar/lib/PBar.o && \
	ar rcs ProgressBar/lib/libPBar.a ProgressBar/lib/PBar.o && \
	$(CXX) -shared -o ProgressBar/lib/PBar.so ProgressBar/lib/PBar.o

# ROOTTools

ROOTToolsLib:
	@$(ECHO) Creating directory for ROOTTools libraries
	mkdir -p ROOTTools/lib

TCanvasPrinter: ROOTTools/src/TCanvasPrinter.cpp ROOTToolsLib
	@$(ECHO) Compiling $@.cpp from ROOTTools module
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

GUIFit: ROOTTools/src/GUIFit.cpp ErrorHandler IOTools ROOTToolsLib
	@$(ECHO) Compiling $@.cpp from ROOTTools module
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) -L./CppTools/lib -l ErrorHandler -lIOTools \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

# SimAnalysis

SimAnalysisLib:
	@$(ECHO) Creating directory for SimAnalysis libraries
	mkdir -p SimAnalysis/lib

EffTreeReader: SimAnalysis/src/EffTreeReader.cpp SimAnalysisLib
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o

EmbTreeReader: SimAnalysis/src/EmbTreeReader.cpp SimAnalysisLib
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o

STrackFun: SimAnalysis/src/STrackFun.cpp SimAnalysisLib
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o SimAnalysis/lib/$@.o && \
	ar rcs SimAnalysis/lib/lib$@.a SimAnalysis/lib/$@.o && \
	$(CXX) -shared -o SimAnalysis/lib/$@.so SimAnalysis/lib/$@.o
	
# SimAnalysis executables

AnalyzeEmbedding: SimAnalysis/src/AnalyzeEmbedding.cpp all_libs bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

AnalyzeHeatMaps: SimAnalysis/src/AnalyzeHeatMaps.cpp all_libs bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

AnalyzeSingleTrack: SimAnalysis/src/AnalyzeSingleTrack.cpp all_libs bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

# other

bin:
	@$(ECHO) Creating directory for binaries
	mkdir bin

clean: 
	@echo Cleaning
	rm -r bin ; \
	rm -r CppTools/lib ; \
	rm -r ROOTTools/lib ; \
	rm -r ProgressBar/lib ; \
	rm -r SimAnalysis/lib

endif
