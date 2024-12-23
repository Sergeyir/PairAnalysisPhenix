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
all_libs: 	 CppToolsLib ROOTToolsLib PBarLib PairAnalysisLib

# CppTools target groups

CppToolsLib: 	 	 ErrorHandler StrTools Time IOTools Box Table
ErrorHandler: 	 	 $(CPP_TOOLS_LIB_DIR)/ErrorHandler.so $(CPP_TOOLS_LIB_DIR)/libErrorHandler.a
StrTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/StrTools.so $(CPP_TOOLS_LIB_DIR)/libStrTools.a
Time: 		  	 	 $(CPP_TOOLS_LIB_DIR)/Time.so $(CPP_TOOLS_LIB_DIR)/libTime.a
IOTools: 	  	 	 $(CPP_TOOLS_LIB_DIR)/IOTools.so $(CPP_TOOLS_LIB_DIR)/libIOTools.a
Box: 			  	 	 $(CPP_TOOLS_LIB_DIR)/Box.so $(CPP_TOOLS_LIB_DIR)/libBox.a
Table: 		  	 	 $(CPP_TOOLS_LIB_DIR)/Table.so $(CPP_TOOLS_LIB_DIR)/libTable.a

# ROOTTools target groups

ROOTToolsLib: 	 	 TCanvasPrinter GUIFit
TCanvasPrinter: 	 $(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.so $(ROOT_TOOLS_LIB_DIR)/libTCanvasPrinter.a
GUIFit: 			 	 $(ROOT_TOOLS_LIB_DIR)/GUIFit.so $(ROOT_TOOLS_LIB_DIR)/libGUIFit.a

# ProgressBar target groups

PBarLib: 		 	 PBar
PBar: 			 	 $(PBAR_LIB_DIR)/PBar.so $(PBAR_LIB_DIR)/libPBar.a

# current repository target groups (Analysis)

PairAnalysisLib: 		 EffTreeReader EmbTreeReader STrackFun PTrackFun IdentFun DataCutsSelector
EffTreeReader:	 	 	 lib/EffTreeReader.so lib/libEffTreeReader.a
EmbTreeReader:	 	 	 lib/EmbTreeReader.so lib/libEmbTreeReader.a
STrackFun:	 	 	 	 lib/STrackFun.so lib/libSTrackFun.a
PTrackFun:	 	 	 	 lib/PTrackFun.so lib/libPTrackFun.a
IdentFun:	 	 	 	 lib/IdentFun.so lib/libIdentFun.a
DataCutsSelector:	 	 lib/DataCutsSelector.so lib/libDataCutsSelector.a

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

$(CPP_TOOLS_LIB_DIR)/IOTools.o: $(CPP_TOOLS_SRC_DIR)/IOTools.cpp \
										  $(CPP_TOOLS_LIB_DIR)/ErrorHandler.o \
										  $(CPP_TOOLS_LIB_DIR)/StrTools.o
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools

$(CPP_TOOLS_LIB_DIR)/Box.o: $(CPP_TOOLS_SRC_DIR)/Box.cpp $(CPP_TOOLS_LIB_DIR)/IOTools.o
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/Table.o: $(CPP_TOOLS_SRC_DIR)/Table.cpp $(CPP_TOOLS_LIB_DIR)/IOTools.o
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) $(CPP_TOOLS_INCLUDE) -o $@ \
	-L $(CPP_TOOLS_LIB_DIR) -lErrorHandler -lStrTools -lIOTools

$(CPP_TOOLS_LIB_DIR)/%.so: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

$(CPP_TOOLS_LIB_DIR)/lib%.a: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX static library $@
	ar rcs $@ $<

# ROOTTools sumbodule targets

$(ROOT_TOOLS_LIB_DIR):
	mkdir -p $@

$(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.o: $(ROOT_TOOLS_SRC_DIR)/TCanvasPrinter.cpp | \
												    $(ROOT_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`

$(ROOT_TOOLS_LIB_DIR)/GUIFit.o: $(ROOT_TOOLS_SRC_DIR)/GUIFit.cpp \
										$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o \
										$(CPP_TOOLS_LIB_DIR)/IOTools.o | \
										$(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) -L$(CPP_TOOLS_LIB_DIR) -l ErrorHandler -lIOTools

$(ROOT_TOOLS_LIB_DIR)/%.so: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

$(ROOT_TOOLS_LIB_DIR)/lib%.a: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Building CXX static library $@
	ar rcs $@ $<

# ProgressBar sumbodule targets

$(PBAR_LIB_DIR):
	mkdir -p $@

$(PBAR_LIB_DIR)/PBar.o: $(PBAR_SRC_DIR)/PBar.cpp | $(PBAR_LIB_DIR)
	@$(ECHO) Building CXX object $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(PBAR_LIB_DIR)/%.so: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

$(PBAR_LIB_DIR)/lib%.a: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Building CXX static library $@
	ar rcs $@ $<

# Current repository targets (Analysis)

lib:
	mkdir -p $@

lib/DataCutsSelector.o: src/DataCutsSelector.cpp | lib
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

lib/%.so: lib/%.o
	@$(ECHO) Building CXX shared library $@
	$(CXX) -shared -o $@ $<

lib/lib%.a: lib/%.o
	@$(ECHO) Building CXX static library $@
	ar rcs $@ $<

AnalyzeHeatMaps: src/AnalyzeHeatMaps.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(PAIR_ANALYSIS_INCLUDE) $(PAIR_ANALYSIS_LIB) 
	
AnalyzeEmbedding: src/AnalyzeEmbedding.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs`
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(PAIR_ANALYSIS_INCLUDE) $(SIM_LIB) \

AnalyzeSingleTrack: src/AnalyzeSingleTrack.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(PAIR_ANALYSIS_INCLUDE) $(SIM_LIB)

AnalyzeResonance: src/AnalyzeResonance.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(PAIR_ANALYSIS_INCLUDE) $(SIM_LIB)

CalibrateSigmalizedResiduals: src/CalibrateSigmalizedResiduals.cpp all_libs bin
	@$(ECHO) Building CXX executable $@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(JSON_INCLUDE) $(JSON_LIB) \
	$(ROOT_INCLUDE) `$(ROOT_CONFIG) --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) $(ROOT_TOOLS_LIB)

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
