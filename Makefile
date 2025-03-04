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
.PHONY: all all_libs exe_targets clean

.SILENT:

# Targets are grouped here

all: all_libs exe_targets
	@echo "All done"

exe_targets: AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack AnalyzeResonance
all_libs: 	 CppToolsLib ROOTToolsLib PBarLib YamlCPP PairAnalysisLib

# CppTools target groups

CppToolsLib: 	 	 ErrorHandler StrTools Time IOTools Box Table
ErrorHandler: 	 	 $(CPP_TOOLS_LIB_DIR)/libErrorHandler.so
StrTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/libStrTools.so
Time: 		  	 	 $(CPP_TOOLS_LIB_DIR)/libTime.so
IOTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/libIOTools.so
Box: 			  	 	 $(CPP_TOOLS_LIB_DIR)/libBox.so
Table: 		  	 	 $(CPP_TOOLS_LIB_DIR)/libTable.so

# ROOTTools target groups

ROOTToolsLib: 	 	 TCanvasPrinter GUIFit
TCanvasPrinter: 	 $(ROOT_TOOLS_LIB_DIR)/libTCanvasPrinter.so
GUIFit: 			 	 $(ROOT_TOOLS_LIB_DIR)/libGUIFit.so

# ProgressBar target groups

PBarLib: 		 	 PBar
PBar: 			 	 $(PBAR_LIB_DIR)/libPBar.so

# yaml-cpp
YamlCPP:			 	 yaml-cpp/build/libyaml-cpp.so

# current repository target groups (Analysis)

PairAnalysisLib: 		 EffTreeReader EmbTreeReader STrackFun PTrackFun \
							 IdentFun DataMethodsSelector InputReader
EffTreeReader:	 	 	 lib/libEffTreeReader.so
EmbTreeReader:	 	 	 lib/libEmbTreeReader.so
STrackFun:	 	 	 	 lib/libSTrackFun.so
PTrackFun:	 	 	 	 lib/libPTrackFun.so
IdentFun:	 	 	 	 lib/libIdentFun.so
DataMethodsSelector:	 lib/libDataMethodsSelector.so
InputReader:	 		 lib/libInputReader.so

# CppTools sumbodule targets

$(CPP_TOOLS_LIB_DIR): 
	mkdir -p $@

$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o: $(CPP_TOOLS_SRC_DIR)/ErrorHandler.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/StrTools.o: $(CPP_TOOLS_SRC_DIR)/StrTools.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/Time.o: $(CPP_TOOLS_SRC_DIR)/Time.cpp | $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(CPP_TOOLS_LIB_DIR)/IOTools.o: $(CPP_TOOLS_SRC_DIR)/IOTools.cpp | ErrorHandler StrTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools

$(CPP_TOOLS_LIB_DIR)/Box.o: $(CPP_TOOLS_SRC_DIR)/Box.cpp | IOTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/Table.o: $(CPP_TOOLS_SRC_DIR)/Table.cpp | IOTools
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/lib%.so: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# ROOTTools sumbodule targets

$(ROOT_TOOLS_LIB_DIR):
	mkdir -p $@

$(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.o: $(ROOT_TOOLS_SRC_DIR)/TCanvasPrinter.cpp | \
												    $(ROOT_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

$(ROOT_TOOLS_LIB_DIR)/GUIFit.o: $(ROOT_TOOLS_SRC_DIR)/GUIFit.cpp | ErrorHandler IOTools \
										  $(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) -L$(CPP_TOOLS_LIB_DIR) -l ErrorHandler -lIOTools

$(ROOT_TOOLS_LIB_DIR)/lib%.so: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# ProgressBar sumbodule targets

$(PBAR_LIB_DIR):
	mkdir -p $@

$(PBAR_LIB_DIR)/PBar.o: $(PBAR_SRC_DIR)/PBar.cpp | $(PBAR_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(PBAR_LIB_DIR)/lib%.so: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

# yaml-cpp

yaml-cpp/build/libyaml-cpp.so:
	@$(ECHO) Building CXX shared library $@
	mkdir -p ./yaml-cpp/build && cmake -S ./yaml-cpp -B ./yaml-cpp/build \
	-DYAML_BUILD_SHARED_LIBS=on && cd ./yaml-cpp/build && make -j4

# Current repository targets (Analysis)

lib:
	mkdir -p $@

lib/InputReader.o: src/InputReader.cpp | lib ErrorHandler YamlCPP
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(JSON_INCLUDE) $(JSON_LIB) \
	$(YAML_INCLUDE) $(YAML_LIB) \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB)

lib/DataMethodsSelector.o: src/DataMethodsSelector.cpp | lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB)

lib/EffTreeReader.o: src/EffTreeReader.cpp | lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

lib/EmbTreeReader.o: src/EmbTreeReader.cpp | lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

lib/STrackFun.o: src/STrackFun.cpp |lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

lib/PTrackFun.o: src/PTrackFun.cpp | lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

lib/IdentFun.o: src/IdentFun.cpp | lib
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ $(PAIR_ANALYSIS_INCLUDE) \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	-L./lib -lSTrackFun

lib/lib%.so: lib/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

AnalyzeHeatMaps: src/AnalyzeHeatMaps.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ALL_INCLUDE) $(ALL_LIB)
	
AnalyzeEmbedding: src/AnalyzeEmbedding.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ALL_INCLUDE) $(ALL_LIB)

AnalyzeSingleTrack: src/AnalyzeSingleTrack.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ALL_INCLUDE) $(ALL_LIB)

AnalyzeResonance: src/AnalyzeResonance.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ALL_INCLUDE) $(ALL_LIB)

# other

bin:
	mkdir $@

clean: 
	@echo Cleaning
	rm -rf bin ; \
	rm -rf CppTools/lib ; \
	rm -rf ROOTTools/lib ; \
	rm -rf ProgressBar/lib ; \
	rm -rf lib

endif
