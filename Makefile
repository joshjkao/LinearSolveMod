CXX = g++-14
CXXFLAGS = -Wall -Wextra -std=c++23

BUILD_DIR = build

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/example example.cpp -lflint

tests: tests.cpp
	$(CXX) $(CXXFLAGS) -I/opt/homebrew/include -L/opt/homebrew/lib -g -o $(BUILD_DIR)/tests tests.cpp -lflint -lgtest -lgtest_main

debug: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o $(BUILD_DIR)/example example.cpp -lflint

run: example
	./$(BUILD_DIR)/example

clean:
	rm -rf $(BUILD_DIR)/*

