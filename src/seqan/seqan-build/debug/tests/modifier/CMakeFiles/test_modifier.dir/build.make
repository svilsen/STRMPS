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
include tests/modifier/CMakeFiles/test_modifier.dir/depend.make

# Include the progress variables for this target.
include tests/modifier/CMakeFiles/test_modifier.dir/progress.make

# Include the compile flags for this target's objects.
include tests/modifier/CMakeFiles/test_modifier.dir/flags.make

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o: tests/modifier/CMakeFiles/test_modifier.dir/flags.make
tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o: /home/svilsen/seqan/seqan-src/tests/modifier/test_modifier.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/modifier && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_modifier.dir/test_modifier.cpp.o -c /home/svilsen/seqan/seqan-src/tests/modifier/test_modifier.cpp

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_modifier.dir/test_modifier.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/modifier && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/modifier/test_modifier.cpp > CMakeFiles/test_modifier.dir/test_modifier.cpp.i

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_modifier.dir/test_modifier.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/modifier && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/modifier/test_modifier.cpp -o CMakeFiles/test_modifier.dir/test_modifier.cpp.s

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.requires:
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.requires

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.provides: tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.requires
	$(MAKE) -f tests/modifier/CMakeFiles/test_modifier.dir/build.make tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.provides.build
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.provides

tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.provides.build: tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o

# Object files for target test_modifier
test_modifier_OBJECTS = \
"CMakeFiles/test_modifier.dir/test_modifier.cpp.o"

# External object files for target test_modifier
test_modifier_EXTERNAL_OBJECTS =

bin/test_modifier: tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o
bin/test_modifier: tests/modifier/CMakeFiles/test_modifier.dir/build.make
bin/test_modifier: tests/modifier/CMakeFiles/test_modifier.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_modifier"
	cd /home/svilsen/seqan/seqan-build/debug/tests/modifier && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_modifier.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/modifier/CMakeFiles/test_modifier.dir/build: bin/test_modifier
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/build

tests/modifier/CMakeFiles/test_modifier.dir/requires: tests/modifier/CMakeFiles/test_modifier.dir/test_modifier.cpp.o.requires
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/requires

tests/modifier/CMakeFiles/test_modifier.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/modifier && $(CMAKE_COMMAND) -P CMakeFiles/test_modifier.dir/cmake_clean.cmake
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/clean

tests/modifier/CMakeFiles/test_modifier.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/modifier /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/modifier /home/svilsen/seqan/seqan-build/debug/tests/modifier/CMakeFiles/test_modifier.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/modifier/CMakeFiles/test_modifier.dir/depend
