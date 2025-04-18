CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++20

BUILD_DIR = build

all: run

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/example example.cpp

gtests: googletests.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/gtests googletests.cpp -lgtest

debug: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o $(BUILD_DIR)/example example.cpp

run: debug
	./$(BUILD_DIR)/example

clean:
	rm -f $(BUILD_DIR)/*

