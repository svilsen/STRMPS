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
include tests/index/CMakeFiles/test_index_finder.dir/depend.make

# Include the progress variables for this target.
include tests/index/CMakeFiles/test_index_finder.dir/progress.make

# Include the compile flags for this target's objects.
include tests/index/CMakeFiles/test_index_finder.dir/flags.make

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o: tests/index/CMakeFiles/test_index_finder.dir/flags.make
tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o: /home/svilsen/seqan/seqan-src/tests/index/test_index_finder.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/index && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o -c /home/svilsen/seqan/seqan-src/tests/index/test_index_finder.cpp

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_index_finder.dir/test_index_finder.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/index && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/index/test_index_finder.cpp > CMakeFiles/test_index_finder.dir/test_index_finder.cpp.i

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_index_finder.dir/test_index_finder.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/index && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/index/test_index_finder.cpp -o CMakeFiles/test_index_finder.dir/test_index_finder.cpp.s

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.requires:
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.requires

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.provides: tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.requires
	$(MAKE) -f tests/index/CMakeFiles/test_index_finder.dir/build.make tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.provides.build
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.provides

tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.provides.build: tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o

# Object files for target test_index_finder
test_index_finder_OBJECTS = \
"CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o"

# External object files for target test_index_finder
test_index_finder_EXTERNAL_OBJECTS =

bin/test_index_finder: tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o
bin/test_index_finder: tests/index/CMakeFiles/test_index_finder.dir/build.make
bin/test_index_finder: tests/index/CMakeFiles/test_index_finder.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_index_finder"
	cd /home/svilsen/seqan/seqan-build/debug/tests/index && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_index_finder.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/index/CMakeFiles/test_index_finder.dir/build: bin/test_index_finder
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/build

tests/index/CMakeFiles/test_index_finder.dir/requires: tests/index/CMakeFiles/test_index_finder.dir/test_index_finder.cpp.o.requires
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/requires

tests/index/CMakeFiles/test_index_finder.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/index && $(CMAKE_COMMAND) -P CMakeFiles/test_index_finder.dir/cmake_clean.cmake
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/clean

tests/index/CMakeFiles/test_index_finder.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/index /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/index /home/svilsen/seqan/seqan-build/debug/tests/index/CMakeFiles/test_index_finder.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/index/CMakeFiles/test_index_finder.dir/depend

