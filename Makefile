CXX = g++-14
CXXFLAGS = -Wall -Wextra -std=c++20

BUILD_DIR = build

all: run

example: example.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/example example.cpp -lflint

gtests: googletests.cpp
	$(CXX) $(CXXFLAGS) -o $(BUILD_DIR)/gtests googletests.cpp -lgtest -lflint

debug: example.cpp
	$(CXX) $(CXXFLAGS) -g -DDEBUG -o $(BUILD_DIR)/example example.cpp -lflint

run: example
	./$(BUILD_DIR)/example

clean:
	rm -f $(BUILD_DIR)/*

