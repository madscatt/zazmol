CXX ?= c++
CC ?= cc
CXXFLAGS ?= -std=c++20 -Wall -Wextra -Wpedantic -O2
CFLAGS ?= -std=c99 -Wall -Wextra -Wpedantic -O2

SASMOL_DIR := ..
BUILD_DIR := $(SASMOL_DIR)/build
LIB := $(BUILD_DIR)/libsasmol_standalone.a
INCLUDES := -I$(SASMOL_DIR)/include
CPP_TARGET := simple_pdb_in_dcd_out_cpp
C_TARGET := simple_pdb_in_dcd_out_c
CPP_SOURCE := simple_pdb_in_dcd_out.cpp
C_SOURCE := simple_pdb_in_dcd_out.c
BRIDGE_SOURCE := sasmol_direct_usage.cpp
CPP_OBJECT := simple_pdb_in_dcd_out_cpp.o
C_OBJECT := simple_pdb_in_dcd_out_c.o
BRIDGE_OBJECT := sasmol_direct_usage.o

INPUT ?= ../../src/python/test_sasmol/data/pdb_common/2AAD.pdb
CPP_OUTPUT ?= 2AAD_from_pdb_cpp.dcd
C_OUTPUT ?= 2AAD_from_pdb_c.dcd

.PHONY: all run run_cpp run_c clean sasmol

all: $(CPP_TARGET) $(C_TARGET)

sasmol:
	cmake -S $(SASMOL_DIR) -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR) --target sasmol_standalone

$(CPP_OBJECT): $(CPP_SOURCE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(CPP_SOURCE) -o $(CPP_OBJECT)

$(C_OBJECT): $(C_SOURCE) sasmol_direct_usage.h
	$(CC) $(CFLAGS) -c $(C_SOURCE) -o $(C_OBJECT)

$(BRIDGE_OBJECT): $(BRIDGE_SOURCE) sasmol_direct_usage.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $(BRIDGE_SOURCE) -o $(BRIDGE_OBJECT)

$(CPP_TARGET): $(CPP_OBJECT) sasmol
	$(CXX) $(CXXFLAGS) $(CPP_OBJECT) $(LIB) -o $(CPP_TARGET)

$(C_TARGET): $(C_OBJECT) $(BRIDGE_OBJECT) sasmol
	$(CXX) $(CXXFLAGS) $(C_OBJECT) $(BRIDGE_OBJECT) $(LIB) -o $(C_TARGET)

run: run_cpp run_c

run_cpp: $(CPP_TARGET)
	./$(CPP_TARGET) $(INPUT) $(CPP_OUTPUT)

run_c: $(C_TARGET)
	./$(C_TARGET) $(INPUT) $(C_OUTPUT)

clean:
	rm -f $(CPP_TARGET) $(C_TARGET) $(CPP_OBJECT) $(C_OBJECT) \
	      $(BRIDGE_OBJECT) $(CPP_OUTPUT) $(C_OUTPUT)
