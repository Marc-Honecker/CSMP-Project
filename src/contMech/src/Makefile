SHELL = /usr/bin/env bash
default_target: seq

COMPILER := g++
CPPFLAGS := -O2 -std=c++11 
#WFLAGS   := -Wno-unused-command-line-argument -Wno-shift-count-overflow -Wno-unused-result
OBJECTS  := auxiliary.o gfmdSheet.o atomicSheet.o interSheet.o contMech.o

# Current date (for backups)
NOW      := $(shell date +%F)
FOLDER   := backup-$(NOW)



# ------ Makefile info ------ #

# This is the reduced Makefile from src15 onward.
# Refer to src14f/Makefile for more legacy options.



# ------ MinGW installation instructions ------ #

# (from http://www.fftw.org/install/windows.html)

# Download latest FFTW3 tar file from http://www.fftw.org/download.html
# Go to the Download folder, unpack, open a MinGW Terminal and execute:
# ./configure
# make
# make install



# ----- Determine correct FFTW paths ----- #

FFTW_PATH_Lnx := /usr/local
FFTW_PATH_Win := C:/MinGW/msys/1.0/local
FFTW_PATH_Vbx := ~/cp2k-master/tools/toolchain/install/fftw-3.3.8
FFTW_PATH_Ssh := ~/bin/fftw

ifneq "$(wildcard ${FFTW_PATH_Win} )" ""
  # If it exists, use the Windows path
  MSG = " Windows FFTW version"
  FFTW_PATH = $(FFTW_PATH_Win)
else ifneq "$(wildcard ${FFTW_PATH_Vbx} )" ""
  # If it exists, use the CP2K toolchain path
  MSG = " CP2K FFTW version"
  FFTW_PATH = $(FFTW_PATH_Vbx)
else ifneq "$(wildcard ${FFTW_PATH_Ssh} )" ""
  # If it exists, use the small cluster path
  MSG = " Cluster FFTW version"
  FFTW_PATH = $(FFTW_PATH_Ssh)
else
  # If none of those exist, use the default path
  MSG = " Linux/Mac default FFTW version"
  FFTW_PATH = $(FFTW_PATH_Lnx)
endif

FFTW_LIB = $(FFTW_PATH)/lib
FFTW_INC = $(FFTW_PATH)/include
LFFTW := -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC)



# ----- contMech compilation rules ----- #

# legacy version
old:
	@rm -f contMech.exe
	@echo $(MSG)
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building contMech.exe ..."
	g++ -O2 contMech.cpp gfmdSheet.cpp \
		-std=c++11 -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC) \
		-o contMech.exe > contMech.e 2>&1

# with optMasses.cpp
opt:
	@rm -f contMech.exe
	@echo $(MSG)
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building contMech.exe ..."
	g++ -O2 contMech.cpp gfmdSheet.cpp optMasses.cpp \
		-std=c++11 -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC) \
		-o contMech.exe > contMech.e 2>&1

dev:
	@echo " 'dev' is deprecated. Use 'make opt' instead."

# with auxiliary.cpp and interSheet.cpp
new:
	@rm -f contMech.exe
	@echo $(MSG)
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " Building contMech.exe ..."
	g++ -O2 contMech.cpp gfmdSheet.cpp atomicSheet.cpp auxiliary.cpp interSheet.cpp \
		-std=c++11 -lfftw3 -L $(FFTW_LIB) -I $(FFTW_INC) \
		-o contMech.exe > contMech.e 2>&1


contMech.o: globals.h contMech.h contMech.cpp
	@echo " compiling contMech.o ..."
	@$(COMPILER) $(CPPFLAGS) $(LFFTW) -c contMech.cpp -o contMech.o > contMech.e 2>&1

gfmdSheet.o: globals.h gfmdSheet.h gfmdSheet.cpp
	@echo " compiling gfmdSheet.o ..."
	@$(COMPILER) $(CPPFLAGS) $(LFFTW) -c gfmdSheet.cpp -o gfmdSheet.o > contMech.e 2>&1

interSheet.o: globals.h interSheet.h interSheet.cpp
	@echo " compiling interSheet.o ..."
	@$(COMPILER) $(CPPFLAGS) $(LFFTW) -c interSheet.cpp -o interSheet.o > contMech.e 2>&1

atomicSheet.o: globals.h atomicSheet.h atomicSheet.cpp
	@echo " compiling atomicSheet.o ..."
	@$(COMPILER) $(CPPFLAGS) $(LFFTW) -c atomicSheet.cpp -o atomicSheet.o > contMech.e 2>&1

auxiliary.o: globals.h auxiliary.cpp
	@echo " compiling auxiliary.o ..."
	@$(COMPILER) $(CPPFLAGS) $(LFFTW) -c auxiliary.cpp -o auxiliary.o > contMech.e 2>&1

seq: $(OBJECTS)
	@echo " ---$(MSG) ---"
	@echo " FFTW Paths: - $(FFTW_LIB)"
	@echo "             - $(FFTW_INC)"
	@echo " compiling contMech.exe ..."
	@$(COMPILER) $(CPPFLAGS) $(OBJECTS) $(LFFTW) -o contMech.exe > contMech.e 2>&1



# ----- cleanup and backup rules ----- #

clean:
	@rm -f *.o
	@rm -f *.exe
	@rm -f *.x
	@rm -f contMech.e


bak:
	@mkdir -p $(FOLDER)
	@cp *.cpp $(FOLDER)
	@cp *.h $(FOLDER)
	@#if [ -f "params.in" ]; then cp -f params.in $(FOLDER); fi
	@#if [ -f "equilPos0.in" ]; then cp -f equilPos0.in $(FOLDER); fi
	@#if [ -f "konfig0E.real" ]; then cp -f konfig0E.real $(FOLDER); fi
	@cp Makefile $(FOLDER)
	@echo " Source code backed up in $(FOLDER)."
