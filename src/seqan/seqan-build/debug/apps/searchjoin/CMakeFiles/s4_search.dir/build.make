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
include apps/searchjoin/CMakeFiles/s4_search.dir/depend.make

# Include the progress variables for this target.
include apps/searchjoin/CMakeFiles/s4_search.dir/progress.make

# Include the compile flags for this target's objects.
include apps/searchjoin/CMakeFiles/s4_search.dir/flags.make

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o: apps/searchjoin/CMakeFiles/s4_search.dir/flags.make
apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o: /home/svilsen/seqan/seqan-src/apps/searchjoin/search.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/searchjoin && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/s4_search.dir/search.cpp.o -c /home/svilsen/seqan/seqan-src/apps/searchjoin/search.cpp

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/s4_search.dir/search.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/searchjoin && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/searchjoin/search.cpp > CMakeFiles/s4_search.dir/search.cpp.i

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/s4_search.dir/search.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/searchjoin && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/searchjoin/search.cpp -o CMakeFiles/s4_search.dir/search.cpp.s

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.requires:
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.requires

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.provides: apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.requires
	$(MAKE) -f apps/searchjoin/CMakeFiles/s4_search.dir/build.make apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.provides.build
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.provides

apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.provides.build: apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o

# Object files for target s4_search
s4_search_OBJECTS = \
"CMakeFiles/s4_search.dir/search.cpp.o"

# External object files for target s4_search
s4_search_EXTERNAL_OBJECTS =

bin/s4_search: apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o
bin/s4_search: apps/searchjoin/CMakeFiles/s4_search.dir/build.make
bin/s4_search: apps/searchjoin/CMakeFiles/s4_search.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/s4_search"
	cd /home/svilsen/seqan/seqan-build/debug/apps/searchjoin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/s4_search.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/searchjoin/CMakeFiles/s4_search.dir/build: bin/s4_search
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/build

apps/searchjoin/CMakeFiles/s4_search.dir/requires: apps/searchjoin/CMakeFiles/s4_search.dir/search.cpp.o.requires
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/requires

apps/searchjoin/CMakeFiles/s4_search.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/searchjoin && $(CMAKE_COMMAND) -P CMakeFiles/s4_search.dir/cmake_clean.cmake
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/clean

apps/searchjoin/CMakeFiles/s4_search.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/searchjoin /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/searchjoin /home/svilsen/seqan/seqan-build/debug/apps/searchjoin/CMakeFiles/s4_search.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/searchjoin/CMakeFiles/s4_search.dir/depend

