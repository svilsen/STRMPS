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
include apps/splazers/CMakeFiles/splazers.dir/depend.make

# Include the progress variables for this target.
include apps/splazers/CMakeFiles/splazers.dir/progress.make

# Include the compile flags for this target's objects.
include apps/splazers/CMakeFiles/splazers.dir/flags.make

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o: apps/splazers/CMakeFiles/splazers.dir/flags.make
apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o: /home/svilsen/seqan/seqan-src/apps/splazers/splazers.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/splazers && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/splazers.dir/splazers.cpp.o -c /home/svilsen/seqan/seqan-src/apps/splazers/splazers.cpp

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/splazers.dir/splazers.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/splazers && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/splazers/splazers.cpp > CMakeFiles/splazers.dir/splazers.cpp.i

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/splazers.dir/splazers.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/splazers && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/splazers/splazers.cpp -o CMakeFiles/splazers.dir/splazers.cpp.s

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.requires:
.PHONY : apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.requires

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.provides: apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.requires
	$(MAKE) -f apps/splazers/CMakeFiles/splazers.dir/build.make apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.provides.build
.PHONY : apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.provides

apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.provides.build: apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o

# Object files for target splazers
splazers_OBJECTS = \
"CMakeFiles/splazers.dir/splazers.cpp.o"

# External object files for target splazers
splazers_EXTERNAL_OBJECTS =

bin/splazers: apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o
bin/splazers: apps/splazers/CMakeFiles/splazers.dir/build.make
bin/splazers: apps/splazers/CMakeFiles/splazers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/splazers"
	cd /home/svilsen/seqan/seqan-build/debug/apps/splazers && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/splazers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/splazers/CMakeFiles/splazers.dir/build: bin/splazers
.PHONY : apps/splazers/CMakeFiles/splazers.dir/build

apps/splazers/CMakeFiles/splazers.dir/requires: apps/splazers/CMakeFiles/splazers.dir/splazers.cpp.o.requires
.PHONY : apps/splazers/CMakeFiles/splazers.dir/requires

apps/splazers/CMakeFiles/splazers.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/splazers && $(CMAKE_COMMAND) -P CMakeFiles/splazers.dir/cmake_clean.cmake
.PHONY : apps/splazers/CMakeFiles/splazers.dir/clean

apps/splazers/CMakeFiles/splazers.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/splazers /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/splazers /home/svilsen/seqan/seqan-build/debug/apps/splazers/CMakeFiles/splazers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/splazers/CMakeFiles/splazers.dir/depend

