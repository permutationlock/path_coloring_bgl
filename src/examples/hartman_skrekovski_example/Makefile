CXX      = g++
CXXFLAGS = -O2 -std=c++1y -Wall -Wno-unused-result
LDFLAGS  =

TARGET = hartman_skrekovski_example
SRCS   = hartman_skrekovski_example.cpp
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)

.PHONY: clean all

all: $(TARGET)

$(TARGET): $(OBJS)
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET)

.cpp.o:
		$(CXX) $(CXXFLAGS) -c $< -o $@

%.depends: %.cpp
		$(CXX) -M $(CXXFLAGS) $< > $@

clean:
		rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)
