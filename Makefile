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

EXE_TARGETS: AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack

ALL_LIBS: $(CPP_TOOLS_LIB_DIR)/ErrorHandler.o $(CPP_TOOLS_LIB_DIR)/ErrorHandler.so \
			 $(CPP_TOOLS_LIB_DIR)/StrTools.o $(CPP_TOOLS_LIB_DIR)/StrTools.so \
			 $(CPP_TOOLS_LIB_DIR)/Time.o $(CPP_TOOLS_LIB_DIR)/Time.so \
			 $(CPP_TOOLS_LIB_DIR)/IOTools.o $(CPP_TOOLS_LIB_DIR)/IOTools.so \
			 $(CPP_TOOLS_LIB_DIR)/Box.o $(CPP_TOOLS_LIB_DIR)/Box.so \
			 $(CPP_TOOLS_LIB_DIR)/Table.o $(CPP_TOOLS_LIB_DIR)/Table.so \
			 $(PBAR_LIB_DIR)/PBar.o $(PBAR_LIB_DIR)/PBar.so \
			 TCanvasPrinter GUIFit EffTreeReader EmbTreeReader STrackFun

.SILENT: $(CPP_TOOLS_LIB_DIR) clean \
			$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o $(CPP_TOOLS_LIB_DIR)/ErrorHandler.so \
			$(CPP_TOOLS_LIB_DIR)/StrTools.o $(CPP_TOOLS_LIB_DIR)/StrTools.so \
			$(CPP_TOOLS_LIB_DIR)/Time.o $(CPP_TOOLS_LIB_DIR)/Time.so \
			$(CPP_TOOLS_LIB_DIR)/IOTools.o $(CPP_TOOLS_LIB_DIR)/IOTools.so \
			$(CPP_TOOLS_LIB_DIR)/Box.o $(CPP_TOOLS_LIB_DIR)/Box.so \
			$(CPP_TOOLS_LIB_DIR)/Table.o $(CPP_TOOLS_LIB_DIR)/Table.so \
			$(PBAR_LIB_DIR)/PBar.o $(PBAR_LIB_DIR)/PBar.so \
			TCanvasPrinter GUIFit EffTreeReader EmbTreeReader STrackFun

all: ALL_LIBS EXE_TARGETS
	@echo "All done"

# CppTools

$(CPP_TOOLS_LIB_DIR): 
	mkdir -p $@

$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o: $(CPP_TOOLS_SRC_DIR)/ErrorHandler.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/StrTools.o: $(CPP_TOOLS_SRC_DIR)/StrTools.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/Time.o: $(CPP_TOOLS_SRC_DIR)/Time.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/IOTools.o: $(CPP_TOOLS_SRC_DIR)/IOTools.cpp \
										  $(CPP_TOOLS_LIB_DIR)/ErrorHandler.o \
										  $(CPP_TOOLS_LIB_DIR)/StrTools.o
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools

$(CPP_TOOLS_LIB_DIR)/Box.o: $(CPP_TOOLS_SRC_DIR)/Box.cpp $(CPP_TOOLS_LIB_DIR)/IOTools.o
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/Table.o: $(CPP_TOOLS_SRC_DIR)/Table.cpp $(CPP_TOOLS_LIB_DIR)/IOTools.o
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/%.so: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(CPP_TOOLS_LIB_DIR)/%.a: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<

# ROOTTools

$(ROOT_TOOLS_LIB_DIR):
	mkdir -p $@

TCanvasPrinter: ROOTTools/src/TCanvasPrinter.cpp ROOTToolsLib
	@$(ECHO) Compiling $@.cpp from ROOTTools submodule
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

GUIFit: ROOTTools/src/GUIFit.cpp ErrorHandler IOTools ROOTToolsLib
	@$(ECHO) Compiling $@.cpp from ROOTTools submodule
	$(CXX) ROOTTools/src/$@.cpp $(CXX_COMMON_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) -L./CppTools/lib -l ErrorHandler -lIOTools \
	-o ROOTTools/lib/$@.o && \
	ar rcs ROOTTools/lib/lib$@.a ROOTTools/lib/$@.o && \
	$(CXX) -shared -o ROOTTools/lib/$@.so ROOTTools/lib/$@.o

# ProgressBar

$(PBAR_LIB_DIR):
	mkdir -p $@

$(PBAR_LIB_DIR)/PBar.o: $(PBAR_SRC_DIR)/PBar.cpp | $(PBAR_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(PBAR_LIB_DIR)/%.so: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(PBAR_LIB_DIR)/%.a: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<


# SimAnalysis

SimAnalysisLib:
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

AnalyzeEmbedding: SimAnalysis/src/AnalyzeEmbedding.cpp ALL_LIBS bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB)

AnalyzeHeatMaps: SimAnalysis/src/AnalyzeHeatMaps.cpp ALL_LIBS bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE)

AnalyzeSingleTrack: SimAnalysis/src/AnalyzeSingleTrack.cpp ALL_LIBS bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE)

AnalyzeResonance: SimAnalysis/src/AnalyzeResonance.cpp ALL_LIBS bin
	@$(ECHO) Compiling $@.cpp from SimAnalysis module
	$(CXX) SimAnalysis/src/$@.cpp -o bin/$@.exe $(CXX_COMMON_EXE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE)

# other

bin:
	mkdir bin

clean: 
	@echo Cleaning
	rm -rf bin/* ; \
	rm -rf CppTools/lib/* ; \
	rm -rf ROOTTools/lib/* ; \
	rm -rf ProgressBar/lib/* ; \
	rm -rf SimAnalysis/lib/*

endif
