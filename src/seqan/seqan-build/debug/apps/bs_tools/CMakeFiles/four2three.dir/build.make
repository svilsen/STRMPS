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
include apps/bs_tools/CMakeFiles/four2three.dir/depend.make

# Include the progress variables for this target.
include apps/bs_tools/CMakeFiles/four2three.dir/progress.make

# Include the compile flags for this target's objects.
include apps/bs_tools/CMakeFiles/four2three.dir/flags.make

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o: apps/bs_tools/CMakeFiles/four2three.dir/flags.make
apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o: /home/svilsen/seqan/seqan-src/apps/bs_tools/four2three.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/bs_tools && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/four2three.dir/four2three.cpp.o -c /home/svilsen/seqan/seqan-src/apps/bs_tools/four2three.cpp

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/four2three.dir/four2three.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/bs_tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/bs_tools/four2three.cpp > CMakeFiles/four2three.dir/four2three.cpp.i

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/four2three.dir/four2three.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/bs_tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/bs_tools/four2three.cpp -o CMakeFiles/four2three.dir/four2three.cpp.s

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.requires:
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.requires

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.provides: apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.requires
	$(MAKE) -f apps/bs_tools/CMakeFiles/four2three.dir/build.make apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.provides.build
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.provides

apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.provides.build: apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o

# Object files for target four2three
four2three_OBJECTS = \
"CMakeFiles/four2three.dir/four2three.cpp.o"

# External object files for target four2three
four2three_EXTERNAL_OBJECTS =

bin/four2three: apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o
bin/four2three: apps/bs_tools/CMakeFiles/four2three.dir/build.make
bin/four2three: /usr/lib/x86_64-linux-gnu/libz.so
bin/four2three: apps/bs_tools/CMakeFiles/four2three.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/four2three"
	cd /home/svilsen/seqan/seqan-build/debug/apps/bs_tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/four2three.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/bs_tools/CMakeFiles/four2three.dir/build: bin/four2three
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/build

apps/bs_tools/CMakeFiles/four2three.dir/requires: apps/bs_tools/CMakeFiles/four2three.dir/four2three.cpp.o.requires
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/requires

apps/bs_tools/CMakeFiles/four2three.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/bs_tools && $(CMAKE_COMMAND) -P CMakeFiles/four2three.dir/cmake_clean.cmake
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/clean

apps/bs_tools/CMakeFiles/four2three.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/bs_tools /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/bs_tools /home/svilsen/seqan/seqan-build/debug/apps/bs_tools/CMakeFiles/four2three.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/bs_tools/CMakeFiles/four2three.dir/depend

