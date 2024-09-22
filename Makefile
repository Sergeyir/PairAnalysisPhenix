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

exe_targets: AnalyzeHeatMaps AnalyzeEmbedding AnalyzeSingleTrack

all_libs: $(CPP_TOOLS_LIB_DIR)/ErrorHandler.o \
			 $(CPP_TOOLS_LIB_DIR)/StrTools.o \
			 $(CPP_TOOLS_LIB_DIR)/Time.o \
			 $(CPP_TOOLS_LIB_DIR)/IOTools.o \
			 $(CPP_TOOLS_LIB_DIR)/Box.o \
			 $(CPP_TOOLS_LIB_DIR)/Table.o \
			 $(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.o \
			 $(ROOT_TOOLS_LIB_DIR)/GUIFit.o \
			 $(PBAR_LIB_DIR)/PBar.o \
			 $(SIM_ANALYSIS_LIB_DIR)/EffTreeReader.o \
			 $(SIM_ANALYSIS_LIB_DIR)/EmbTreeReader.o \
			 $(SIM_ANALYSIS_LIB_DIR)/STrackFun.o \
			 $(SIM_ANALYSIS_LIB_DIR)/PTrackFun.o \
			 $(SIM_ANALYSIS_LIB_DIR)/IdentFun.o \
			 $(ANALYSIS_LIB_DIR)/DeadAreasCuts.o \
			 $(CPP_TOOLS_LIB_DIR)/ErrorHandler.so \
			 $(CPP_TOOLS_LIB_DIR)/StrTools.so \
			 $(CPP_TOOLS_LIB_DIR)/Time.so \
			 $(CPP_TOOLS_LIB_DIR)/IOTools.so \
			 $(CPP_TOOLS_LIB_DIR)/Box.so \
			 $(CPP_TOOLS_LIB_DIR)/Table.so \
			 $(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.so \
			 $(ROOT_TOOLS_LIB_DIR)/GUIFit.so \
			 $(PBAR_LIB_DIR)/PBar.so \
			 $(SIM_ANALYSIS_LIB_DIR)/EffTreeReader.so \
			 $(SIM_ANALYSIS_LIB_DIR)/EmbTreeReader.so \
			 $(SIM_ANALYSIS_LIB_DIR)/STrackFun.so \
			 $(SIM_ANALYSIS_LIB_DIR)/PTrackFun.so \
			 $(SIM_ANALYSIS_LIB_DIR)/IdentFun.so \
			 $(ANALYSIS_LIB_DIR)/DeadAreasCuts.so \
			 $(CPP_TOOLS_LIB_DIR)/libErrorHandler.a \
			 $(CPP_TOOLS_LIB_DIR)/libStrTools.a \
			 $(CPP_TOOLS_LIB_DIR)/libTime.a \
			 $(CPP_TOOLS_LIB_DIR)/libIOTools.a \
			 $(CPP_TOOLS_LIB_DIR)/libBox.a \
			 $(CPP_TOOLS_LIB_DIR)/libTable.a \
			 $(ROOT_TOOLS_LIB_DIR)/libTCanvasPrinter.a \
			 $(ROOT_TOOLS_LIB_DIR)/libGUIFit.a \
			 $(PBAR_LIB_DIR)/libPBar.a \
			 $(SIM_ANALYSIS_LIB_DIR)/libEffTreeReader.a \
			 $(SIM_ANALYSIS_LIB_DIR)/libEmbTreeReader.a \
			 $(SIM_ANALYSIS_LIB_DIR)/libSTrackFun.a \
			 $(SIM_ANALYSIS_LIB_DIR)/libPTrackFun.a \
			 $(SIM_ANALYSIS_LIB_DIR)/libIdentFun.a \
			 $(ANALYSIS_LIB_DIR)/libDeadAreasCuts.a

.SILENT:

all: all_libs exe_targets
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

$(CPP_TOOLS_LIB_DIR)/lib%.a: $(CPP_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<

# ROOTTools

$(ROOT_TOOLS_LIB_DIR):
	mkdir -p $@

$(ROOT_TOOLS_LIB_DIR)/TCanvasPrinter.o: $(ROOT_TOOLS_SRC_DIR)/TCanvasPrinter.cpp | \
												    $(ROOT_TOOLS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

$(ROOT_TOOLS_LIB_DIR)/GUIFit.o: $(ROOT_TOOLS_SRC_DIR)/GUIFit.cpp \
										$(CPP_TOOLS_LIB_DIR)/ErrorHandler.o \
										$(CPP_TOOLS_LIB_DIR)/IOTools.o | \
										$(CPP_TOOLS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) -L$(CPP_TOOLS_LIB_DIR) -l ErrorHandler -lIOTools

$(ROOT_TOOLS_LIB_DIR)/%.so: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(ROOT_TOOLS_LIB_DIR)/lib%.a: $(ROOT_TOOLS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<

# ProgressBar

$(PBAR_LIB_DIR):
	mkdir -p $@

$(PBAR_LIB_DIR)/PBar.o: $(PBAR_SRC_DIR)/PBar.cpp | $(PBAR_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(PBAR_LIB_DIR)/%.so: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(PBAR_LIB_DIR)/lib%.a: $(PBAR_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<


# SimAnalysis

$(SIM_ANALYSIS_LIB_DIR):
	mkdir -p $@

$(SIM_ANALYSIS_LIB_DIR)/EffTreeReader.o: $(SIM_ANALYSIS_SRC_DIR)/EffTreeReader.cpp | \
													  $(SIM_ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

$(SIM_ANALYSIS_LIB_DIR)/EmbTreeReader.o: $(SIM_ANALYSIS_SRC_DIR)/EmbTreeReader.cpp | \
													  $(SIM_ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

$(SIM_ANALYSIS_LIB_DIR)/STrackFun.o: $(SIM_ANALYSIS_SRC_DIR)/STrackFun.cpp | \
											    $(SIM_ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

$(SIM_ANALYSIS_LIB_DIR)/PTrackFun.o: $(SIM_ANALYSIS_SRC_DIR)/PTrackFun.cpp | \
											    $(SIM_ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

$(SIM_ANALYSIS_LIB_DIR)/IdentFun.o: $(SIM_ANALYSIS_SRC_DIR)/IdentFun.cpp | \
											   $(SIM_ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@ $(ANALYSIS_INCLUDE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(SIM_ANALYSIS_INCLUDE) -L./$(SIM_ANALYSIS_LIB_DIR) -lSTrackFun \
	$(ANALYSIS_INCLUDE) -L./$(ANALYSIS_LIB_DIR) -lDeadAreasCuts

$(SIM_ANALYSIS_LIB_DIR)/%.so: $(SIM_ANALYSIS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(SIM_ANALYSIS_LIB_DIR)/lib%.a: $(SIM_ANALYSIS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<
	
# SimAnalysis executables

AnalyzeEmbedding: SimAnalysis/src/AnalyzeEmbedding.cpp all_libs bin
	@$(ECHO) Compiling $< into bin/$@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`

AnalyzeHeatMaps: SimAnalysis/src/AnalyzeHeatMaps.cpp all_libs bin
	@$(ECHO) Compiling $< into bin/$@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE)

AnalyzeSingleTrack: SimAnalysis/src/AnalyzeSingleTrack.cpp all_libs bin
	@$(ECHO) Compiling $< into bin/$@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE)

AnalyzeResonance: SimAnalysis/src/AnalyzeResonance.cpp all_libs bin
	@$(ECHO) Compiling $< into bin/$@
	$(CXX) $< $(CXX_COMMON_EXE) -o bin/$@ \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs` \
	$(CPP_TOOLS_INCLUDE) $(CPP_TOOLS_LIB) \
	$(PBAR_INCLUDE) $(PBAR_LIB) \
	$(ROOT_TOOLS_INCLUDE) \
	$(SIM_ANALYSIS_INCLUDE) $(SIM_ANALYSIS_LIB) \
	$(ANALYSIS_INCLUDE) $(ANALYSIS_LIB)

# Analysis

$(ANALYSIS_LIB_DIR):
	mkdir -p $@

$(ANALYSIS_LIB_DIR)/DeadAreasCuts.o: $(ANALYSIS_SRC_DIR)/DeadAreasCuts.cpp | $(ANALYSIS_LIB_DIR)
	@$(ECHO) Compiling $< into $@
	$(CXX) $< $(CXX_COMMON_LIB) -o $@

$(ANALYSIS_LIB_DIR)/%.so: $(ANALYSIS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	$(CXX) -shared -o $@ $<

$(ANALYSIS_LIB_DIR)/lib%.a: $(ANALYSIS_LIB_DIR)/%.o
	@$(ECHO) Creating $@ from $<
	ar rcs $@ $<

# other

bin:
	mkdir $@

clean: 
	@echo Cleaning
	rm -rf bin ; \
	rm -rf CppTools/lib ; \
	rm -rf ROOTTools/lib ; \
	rm -rf ProgressBar/lib ; \
	rm -rf SimAnalysis/lib

endif
