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
include apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/depend.make

# Include the progress variables for this target.
include apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/progress.make

# Include the compile flags for this target's objects.
include apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/flags.make

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/flags.make
apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o: /home/svilsen/seqan/seqan-src/apps/pair_align/lib/pair_align_global.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o -c /home/svilsen/seqan/seqan-src/apps/pair_align/lib/pair_align_global.cpp

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/pair_align/lib/pair_align_global.cpp > CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.i

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/pair_align/lib/pair_align_global.cpp -o CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.s

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.requires:
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.requires

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.provides: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.requires
	$(MAKE) -f apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/build.make apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.provides.build
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.provides

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.provides.build: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o

# Object files for target pair_align_global_1011
pair_align_global_1011_OBJECTS = \
"CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o"

# External object files for target pair_align_global_1011
pair_align_global_1011_EXTERNAL_OBJECTS =

apps/pair_align/lib/libpair_align_global_1011.a: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o
apps/pair_align/lib/libpair_align_global_1011.a: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/build.make
apps/pair_align/lib/libpair_align_global_1011.a: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libpair_align_global_1011.a"
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && $(CMAKE_COMMAND) -P CMakeFiles/pair_align_global_1011.dir/cmake_clean_target.cmake
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pair_align_global_1011.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/build: apps/pair_align/lib/libpair_align_global_1011.a
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/build

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/requires: apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/pair_align_global.cpp.o.requires
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/requires

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib && $(CMAKE_COMMAND) -P CMakeFiles/pair_align_global_1011.dir/cmake_clean.cmake
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/clean

apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/pair_align/lib /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib /home/svilsen/seqan/seqan-build/debug/apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/pair_align/lib/CMakeFiles/pair_align_global_1011.dir/depend

