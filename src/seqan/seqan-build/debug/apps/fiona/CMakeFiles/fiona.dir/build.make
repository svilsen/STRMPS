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
include apps/fiona/CMakeFiles/fiona.dir/depend.make

# Include the progress variables for this target.
include apps/fiona/CMakeFiles/fiona.dir/progress.make

# Include the compile flags for this target's objects.
include apps/fiona/CMakeFiles/fiona.dir/flags.make

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o: apps/fiona/CMakeFiles/fiona.dir/flags.make
apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o: /home/svilsen/seqan/seqan-src/apps/fiona/fiona.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fiona && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fiona.dir/fiona.cpp.o -c /home/svilsen/seqan/seqan-src/apps/fiona/fiona.cpp

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fiona.dir/fiona.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fiona && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/fiona/fiona.cpp > CMakeFiles/fiona.dir/fiona.cpp.i

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fiona.dir/fiona.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fiona && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/fiona/fiona.cpp -o CMakeFiles/fiona.dir/fiona.cpp.s

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.requires:
.PHONY : apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.requires

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.provides: apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.requires
	$(MAKE) -f apps/fiona/CMakeFiles/fiona.dir/build.make apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.provides.build
.PHONY : apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.provides

apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.provides.build: apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o

# Object files for target fiona
fiona_OBJECTS = \
"CMakeFiles/fiona.dir/fiona.cpp.o"

# External object files for target fiona
fiona_EXTERNAL_OBJECTS =

bin/fiona: apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o
bin/fiona: apps/fiona/CMakeFiles/fiona.dir/build.make
bin/fiona: /usr/lib/x86_64-linux-gnu/libz.so
bin/fiona: apps/fiona/CMakeFiles/fiona.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/fiona"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fiona && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fiona.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/fiona/CMakeFiles/fiona.dir/build: bin/fiona
.PHONY : apps/fiona/CMakeFiles/fiona.dir/build

apps/fiona/CMakeFiles/fiona.dir/requires: apps/fiona/CMakeFiles/fiona.dir/fiona.cpp.o.requires
.PHONY : apps/fiona/CMakeFiles/fiona.dir/requires

apps/fiona/CMakeFiles/fiona.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/fiona && $(CMAKE_COMMAND) -P CMakeFiles/fiona.dir/cmake_clean.cmake
.PHONY : apps/fiona/CMakeFiles/fiona.dir/clean

apps/fiona/CMakeFiles/fiona.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/fiona /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/fiona /home/svilsen/seqan/seqan-build/debug/apps/fiona/CMakeFiles/fiona.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/fiona/CMakeFiles/fiona.dir/depend

