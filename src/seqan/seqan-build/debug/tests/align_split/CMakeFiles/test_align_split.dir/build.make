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
include tests/align_split/CMakeFiles/test_align_split.dir/depend.make

# Include the progress variables for this target.
include tests/align_split/CMakeFiles/test_align_split.dir/progress.make

# Include the compile flags for this target's objects.
include tests/align_split/CMakeFiles/test_align_split.dir/flags.make

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o: tests/align_split/CMakeFiles/test_align_split.dir/flags.make
tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o: /home/svilsen/seqan/seqan-src/tests/align_split/test_align_split.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/align_split && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_align_split.dir/test_align_split.cpp.o -c /home/svilsen/seqan/seqan-src/tests/align_split/test_align_split.cpp

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_align_split.dir/test_align_split.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/align_split && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/align_split/test_align_split.cpp > CMakeFiles/test_align_split.dir/test_align_split.cpp.i

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_align_split.dir/test_align_split.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/align_split && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/align_split/test_align_split.cpp -o CMakeFiles/test_align_split.dir/test_align_split.cpp.s

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.requires:
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.requires

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.provides: tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.requires
	$(MAKE) -f tests/align_split/CMakeFiles/test_align_split.dir/build.make tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.provides.build
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.provides

tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.provides.build: tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o

# Object files for target test_align_split
test_align_split_OBJECTS = \
"CMakeFiles/test_align_split.dir/test_align_split.cpp.o"

# External object files for target test_align_split
test_align_split_EXTERNAL_OBJECTS =

bin/test_align_split: tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o
bin/test_align_split: tests/align_split/CMakeFiles/test_align_split.dir/build.make
bin/test_align_split: tests/align_split/CMakeFiles/test_align_split.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_align_split"
	cd /home/svilsen/seqan/seqan-build/debug/tests/align_split && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_align_split.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/align_split/CMakeFiles/test_align_split.dir/build: bin/test_align_split
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/build

tests/align_split/CMakeFiles/test_align_split.dir/requires: tests/align_split/CMakeFiles/test_align_split.dir/test_align_split.cpp.o.requires
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/requires

tests/align_split/CMakeFiles/test_align_split.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/align_split && $(CMAKE_COMMAND) -P CMakeFiles/test_align_split.dir/cmake_clean.cmake
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/clean

tests/align_split/CMakeFiles/test_align_split.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/align_split /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/align_split /home/svilsen/seqan/seqan-build/debug/tests/align_split/CMakeFiles/test_align_split.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/align_split/CMakeFiles/test_align_split.dir/depend

