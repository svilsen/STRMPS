# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/svilsen/seqan/seqan-src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/svilsen/seqan/seqan-build/debug

# Include any dependencies generated for this target.
include tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/depend.make

# Include the progress variables for this target.
include tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/progress.make

# Include the compile flags for this target's objects.
include tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/flags.make

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/flags.make
tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o: /home/svilsen/seqan/seqan-src/tests/graph_algorithms/test_graph_algorithms.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o -c /home/svilsen/seqan/seqan-src/tests/graph_algorithms/test_graph_algorithms.cpp

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/graph_algorithms/test_graph_algorithms.cpp > CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.i

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/graph_algorithms/test_graph_algorithms.cpp -o CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.s

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.requires:
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.requires

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.provides: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.requires
	$(MAKE) -f tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/build.make tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.provides.build
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.provides

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.provides.build: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o

# Object files for target test_graph_algorithms
test_graph_algorithms_OBJECTS = \
"CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o"

# External object files for target test_graph_algorithms
test_graph_algorithms_EXTERNAL_OBJECTS =

bin/test_graph_algorithms: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o
bin/test_graph_algorithms: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/build.make
bin/test_graph_algorithms: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_graph_algorithms"
	cd /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_graph_algorithms.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/build: bin/test_graph_algorithms
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/build

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/requires: tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/test_graph_algorithms.cpp.o.requires
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/requires

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms && $(CMAKE_COMMAND) -P CMakeFiles/test_graph_algorithms.dir/cmake_clean.cmake
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/clean

tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/graph_algorithms /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms /home/svilsen/seqan/seqan-build/debug/tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/graph_algorithms/CMakeFiles/test_graph_algorithms.dir/depend

