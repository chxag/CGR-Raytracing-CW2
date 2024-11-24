# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -Wall -Wextra -I/opt/homebrew/include

# Include directories for headers (if you have headers in 'include' folder)
INCLUDES = -Iinclude

# Source files
SRCS = raytracer.cpp tools.cpp sphere.cpp ppmWriter.cpp triangle.cpp cylinder.cpp blinn_phong_shader.cpp binary_shader.cpp vector_utils.cpp shadow.cpp tone_mapping.cpp texture.cpp

# Header files (add header files if needed for dependencies)
HDRS = sphere.h tools.h triangle.h cylinder.h blinn_phong_shader.h binary_shader.h vector_utils.h shadow.h tone_mapping.h texture.h

# Target executable
TARGET = raytracer

# Object files (replace .cpp with .o)
OBJS = $(SRCS:.cpp=.o)

# Default target (build the program)
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean