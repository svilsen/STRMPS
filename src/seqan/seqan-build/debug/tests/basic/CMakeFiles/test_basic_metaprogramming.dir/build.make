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
include tests/basic/CMakeFiles/test_basic_metaprogramming.dir/depend.make

# Include the progress variables for this target.
include tests/basic/CMakeFiles/test_basic_metaprogramming.dir/progress.make

# Include the compile flags for this target's objects.
include tests/basic/CMakeFiles/test_basic_metaprogramming.dir/flags.make

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/flags.make
tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o: /home/svilsen/seqan/seqan-src/tests/basic/test_basic_metaprogramming.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/basic && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o -c /home/svilsen/seqan/seqan-src/tests/basic/test_basic_metaprogramming.cpp

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/basic && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/basic/test_basic_metaprogramming.cpp > CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.i

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/basic && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/basic/test_basic_metaprogramming.cpp -o CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.s

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.requires:
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.requires

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.provides: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.requires
	$(MAKE) -f tests/basic/CMakeFiles/test_basic_metaprogramming.dir/build.make tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.provides.build
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.provides

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.provides.build: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o

# Object files for target test_basic_metaprogramming
test_basic_metaprogramming_OBJECTS = \
"CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o"

# External object files for target test_basic_metaprogramming
test_basic_metaprogramming_EXTERNAL_OBJECTS =

bin/test_basic_metaprogramming: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o
bin/test_basic_metaprogramming: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/build.make
bin/test_basic_metaprogramming: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_basic_metaprogramming"
	cd /home/svilsen/seqan/seqan-build/debug/tests/basic && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_basic_metaprogramming.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/basic/CMakeFiles/test_basic_metaprogramming.dir/build: bin/test_basic_metaprogramming
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/build

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/requires: tests/basic/CMakeFiles/test_basic_metaprogramming.dir/test_basic_metaprogramming.cpp.o.requires
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/requires

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/basic && $(CMAKE_COMMAND) -P CMakeFiles/test_basic_metaprogramming.dir/cmake_clean.cmake
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/clean

tests/basic/CMakeFiles/test_basic_metaprogramming.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/basic /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/basic /home/svilsen/seqan/seqan-build/debug/tests/basic/CMakeFiles/test_basic_metaprogramming.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/basic/CMakeFiles/test_basic_metaprogramming.dir/depend

