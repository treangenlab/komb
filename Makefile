CXX := g++
CXXFLAGS := -O3 -fopenmp -std=c++11 -Wl,-rpath=$(shell pwd)/lib/
LDFLAGS = -Llib/ -ligraph 
BUILD := ./build
OBJ_DIR := $(BUILD)/objects
APP_DIR := $(BUILD)/apps
TARGET := komb
INCLUDE_IGRAPH := -Iinclude/igraph/
INLCUDE_BOOST := -Iinclude/boost
INLCUDE_MISC := -Iinclude/
SRC := $(wildcard src/*.cpp)

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDE_IGRAPH) $(INCLUDE_BOOST) $(INLCUDE_MISC)  -o $@  

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE_IGRAPH) $(INLCUDE_BOOST) $(INCLUDE_MISC) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LDFLAGS)

.PHONY: all build clean

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

