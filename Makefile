CXX = g++
CXXFLAGS = -std=c++14 -O2 -Wall -I/usr/include/eigen3

SRC_DIR = src
INCLUDE_DIR = include
SCRIPTS_DIR = scripts
BUILD_DIR = build

# Source files
SOURCES = $(SRC_DIR)/InteriorPointLP.cpp

# Test executables
TEST_DEBUG = $(BUILD_DIR)/test_debug
TEST3_MAIN = $(BUILD_DIR)/test3_main

all: directories $(TEST_DEBUG) $(TEST3_MAIN)

directories:
	mkdir -p $(BUILD_DIR)

$(TEST_DEBUG): $(SCRIPTS_DIR)/test_debug.cpp $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $^ -o $@

$(TEST3_MAIN): $(SCRIPTS_DIR)/test3_main.cpp $(SOURCES)
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) $^ -o $@

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all directories clean
