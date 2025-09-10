CXX = g++-14
CXXFLAGS = -Wall -Wextra -Wpedantic -std=c++23

BUILD_DIR = build

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/example example.cpp

gtests: googletests.cpp
	$(CXX) $(CXXFLAGS) -I/opt/homebrew/include -L/opt/homebrew/lib -g -o $(BUILD_DIR)/gtests googletests.cpp -lgtest -lgtest_main

debug: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o $(BUILD_DIR)/example example.cpp

flint: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -DFLINT -I/opt/homebrew/include -L/opt/homebrew/lib -o $(BUILD_DIR)/example example.cpp -lflint

run: example
	./$(BUILD_DIR)/example

clean:
	rm -rf $(BUILD_DIR)/*

