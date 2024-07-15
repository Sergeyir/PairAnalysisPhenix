# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

all: ErrorHandler StrTools Time

ErrorHandler: CppTools/src/ErrorHandler.cpp lib
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o
	
StrTools: CppTools/src/StrTools.cpp lib
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o

Time: CppTools/src/Time.cpp lib
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o

IOTools: CppTools/src/IOTools.cpp ErrorHandler StrTools
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	$(CPP_TOOLS_INCLUDE) \
	$(ERROR_HANDLER_LIB) $(STR_TOOLS_LIB) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o

Box: CppTools/src/Box.cpp IOTools
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	$(CPP_TOOLS_INCLUDE) \
	$(ERROR_HANDLER_LIB) $(STR_TOOLS_LIB) $(IOTOOLS_LIB) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o

Table: CppTools/src/Table.cpp IOTools
	$(CXX) CppTools/src/$@.cpp $(CXX_FLAGS) \
	$(CPP_TOOLS_INCLUDE) \
	$(ERROR_HANDLER_LIB) $(STR_TOOLS_LIB) $(IOTOOLS_LIB) \
	-o lib/$@.o && \
	ar rcs lib/lib$@.a lib/$@.o && \
	$(CXX) -shared -o lib/$@.so lib/$@.o

lib: 
	mkdir lib

clean: 
	rm -rf lib/*
