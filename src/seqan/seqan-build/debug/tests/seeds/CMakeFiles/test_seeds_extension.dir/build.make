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
include tests/seeds/CMakeFiles/test_seeds_extension.dir/depend.make

# Include the progress variables for this target.
include tests/seeds/CMakeFiles/test_seeds_extension.dir/progress.make

# Include the compile flags for this target's objects.
include tests/seeds/CMakeFiles/test_seeds_extension.dir/flags.make

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o: tests/seeds/CMakeFiles/test_seeds_extension.dir/flags.make
tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o: /home/svilsen/seqan/seqan-src/tests/seeds/test_seeds_extension.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/tests/seeds && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o -c /home/svilsen/seqan/seqan-src/tests/seeds/test_seeds_extension.cpp

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/tests/seeds && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/tests/seeds/test_seeds_extension.cpp > CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.i

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/tests/seeds && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/tests/seeds/test_seeds_extension.cpp -o CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.s

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.requires:
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.requires

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.provides: tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.requires
	$(MAKE) -f tests/seeds/CMakeFiles/test_seeds_extension.dir/build.make tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.provides.build
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.provides

tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.provides.build: tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o

# Object files for target test_seeds_extension
test_seeds_extension_OBJECTS = \
"CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o"

# External object files for target test_seeds_extension
test_seeds_extension_EXTERNAL_OBJECTS =

bin/test_seeds_extension: tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o
bin/test_seeds_extension: tests/seeds/CMakeFiles/test_seeds_extension.dir/build.make
bin/test_seeds_extension: tests/seeds/CMakeFiles/test_seeds_extension.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_seeds_extension"
	cd /home/svilsen/seqan/seqan-build/debug/tests/seeds && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_seeds_extension.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/seeds/CMakeFiles/test_seeds_extension.dir/build: bin/test_seeds_extension
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/build

tests/seeds/CMakeFiles/test_seeds_extension.dir/requires: tests/seeds/CMakeFiles/test_seeds_extension.dir/test_seeds_extension.cpp.o.requires
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/requires

tests/seeds/CMakeFiles/test_seeds_extension.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/tests/seeds && $(CMAKE_COMMAND) -P CMakeFiles/test_seeds_extension.dir/cmake_clean.cmake
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/clean

tests/seeds/CMakeFiles/test_seeds_extension.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/tests/seeds /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/tests/seeds /home/svilsen/seqan/seqan-build/debug/tests/seeds/CMakeFiles/test_seeds_extension.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/seeds/CMakeFiles/test_seeds_extension.dir/depend

