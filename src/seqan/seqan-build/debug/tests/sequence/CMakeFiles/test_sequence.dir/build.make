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
include tests/sequence/CMakeFiles/test_sequence.dir/depend.make

# Include the progress variables for this target.
include tests/sequence/CMakeFiles/test_sequence.dir/progress.make

# Include the compile flags for this target's objects.
include tests/sequence/CMakeFiles/test_sequence.dir/flags.make

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o: tests/sequence/CMakeFiles/test_sequence.dir/flags.make
tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o: /home/svilsen/seqan/seqan-src/tests/sequence/test_sequence.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_sequence.dir/test_sequence.cpp.o -c /home/svilsen/seqan/seqan-src/tests/sequence/test_sequence.cpp

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_sequence.dir/test_sequence.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/sequence/test_sequence.cpp > CMakeFiles/test_sequence.dir/test_sequence.cpp.i

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_sequence.dir/test_sequence.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/sequence/test_sequence.cpp -o CMakeFiles/test_sequence.dir/test_sequence.cpp.s

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.requires:
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.requires

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.provides: tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.requires
	$(MAKE) -f tests/sequence/CMakeFiles/test_sequence.dir/build.make tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.provides.build
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.provides

tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.provides.build: tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o

# Object files for target test_sequence
test_sequence_OBJECTS = \
"CMakeFiles/test_sequence.dir/test_sequence.cpp.o"

# External object files for target test_sequence
test_sequence_EXTERNAL_OBJECTS =

bin/test_sequence: tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o
bin/test_sequence: tests/sequence/CMakeFiles/test_sequence.dir/build.make
bin/test_sequence: tests/sequence/CMakeFiles/test_sequence.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_sequence"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sequence.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/sequence/CMakeFiles/test_sequence.dir/build: bin/test_sequence
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/build

tests/sequence/CMakeFiles/test_sequence.dir/requires: tests/sequence/CMakeFiles/test_sequence.dir/test_sequence.cpp.o.requires
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/requires

tests/sequence/CMakeFiles/test_sequence.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence && $(CMAKE_COMMAND) -P CMakeFiles/test_sequence.dir/cmake_clean.cmake
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/clean

tests/sequence/CMakeFiles/test_sequence.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/sequence /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/sequence /home/svilsen/seqan/seqan-build/debug/tests/sequence/CMakeFiles/test_sequence.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/sequence/CMakeFiles/test_sequence.dir/depend

