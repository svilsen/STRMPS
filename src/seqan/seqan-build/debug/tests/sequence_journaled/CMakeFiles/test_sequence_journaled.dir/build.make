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
include tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/depend.make

# Include the progress variables for this target.
include tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/progress.make

# Include the compile flags for this target's objects.
include tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/flags.make

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/flags.make
tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o: /home/svilsen/seqan/seqan-src/tests/sequence_journaled/test_sequence_journaled.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o -c /home/svilsen/seqan/seqan-src/tests/sequence_journaled/test_sequence_journaled.cpp

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/sequence_journaled/test_sequence_journaled.cpp > CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.i

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/sequence_journaled/test_sequence_journaled.cpp -o CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.s

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.requires:
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.requires

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.provides: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.requires
	$(MAKE) -f tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/build.make tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.provides.build
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.provides

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.provides.build: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o

# Object files for target test_sequence_journaled
test_sequence_journaled_OBJECTS = \
"CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o"

# External object files for target test_sequence_journaled
test_sequence_journaled_EXTERNAL_OBJECTS =

bin/test_sequence_journaled: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o
bin/test_sequence_journaled: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/build.make
bin/test_sequence_journaled: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_sequence_journaled"
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sequence_journaled.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/build: bin/test_sequence_journaled
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/build

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/requires: tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/test_sequence_journaled.cpp.o.requires
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/requires

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled && $(CMAKE_COMMAND) -P CMakeFiles/test_sequence_journaled.dir/cmake_clean.cmake
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/clean

tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/sequence_journaled /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled /home/svilsen/seqan/seqan-build/debug/tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/sequence_journaled/CMakeFiles/test_sequence_journaled.dir/depend

