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
include tests/parse_lm/CMakeFiles/test_parse_lm.dir/depend.make

# Include the progress variables for this target.
include tests/parse_lm/CMakeFiles/test_parse_lm.dir/progress.make

# Include the compile flags for this target's objects.
include tests/parse_lm/CMakeFiles/test_parse_lm.dir/flags.make

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o: tests/parse_lm/CMakeFiles/test_parse_lm.dir/flags.make
tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o: /home/svilsen/seqan/seqan-src/tests/parse_lm/test_parse_lm.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/parse_lm && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o -c /home/svilsen/seqan/seqan-src/tests/parse_lm/test_parse_lm.cpp

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/parse_lm && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/parse_lm/test_parse_lm.cpp > CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.i

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/parse_lm && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/parse_lm/test_parse_lm.cpp -o CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.s

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.requires:
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.requires

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.provides: tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.requires
	$(MAKE) -f tests/parse_lm/CMakeFiles/test_parse_lm.dir/build.make tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.provides.build
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.provides

tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.provides.build: tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o

# Object files for target test_parse_lm
test_parse_lm_OBJECTS = \
"CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o"

# External object files for target test_parse_lm
test_parse_lm_EXTERNAL_OBJECTS =

bin/test_parse_lm: tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o
bin/test_parse_lm: tests/parse_lm/CMakeFiles/test_parse_lm.dir/build.make
bin/test_parse_lm: tests/parse_lm/CMakeFiles/test_parse_lm.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_parse_lm"
	cd /home/svilsen/seqan/seqan-build/debug/tests/parse_lm && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_parse_lm.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/parse_lm/CMakeFiles/test_parse_lm.dir/build: bin/test_parse_lm
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/build

tests/parse_lm/CMakeFiles/test_parse_lm.dir/requires: tests/parse_lm/CMakeFiles/test_parse_lm.dir/test_parse_lm.cpp.o.requires
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/requires

tests/parse_lm/CMakeFiles/test_parse_lm.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/parse_lm && $(CMAKE_COMMAND) -P CMakeFiles/test_parse_lm.dir/cmake_clean.cmake
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/clean

tests/parse_lm/CMakeFiles/test_parse_lm.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/parse_lm /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/parse_lm /home/svilsen/seqan/seqan-build/debug/tests/parse_lm/CMakeFiles/test_parse_lm.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/parse_lm/CMakeFiles/test_parse_lm.dir/depend

