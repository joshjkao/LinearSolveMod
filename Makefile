CXX = g++-14
CXXFLAGS = -Wall -Wextra -std=c++23

BUILD_DIR = build

all: run

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/example example.cpp

gtests: googletests.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/gtests googletests.cpp -lgtest

debug: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o $(BUILD_DIR)/example example.cpp

flint: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -DFLINT -I/opt/homebrew/include -L/opt/homebrew/lib -o $(BUILD_DIR)/example example.cpp -lflint

run: example
	./$(BUILD_DIR)/example

clean:
	rm -rf $(BUILD_DIR)/*

