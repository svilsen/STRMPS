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
include apps/mason2/CMakeFiles/mason_materializer.dir/depend.make

# Include the progress variables for this target.
include apps/mason2/CMakeFiles/mason_materializer.dir/progress.make

# Include the compile flags for this target's objects.
include apps/mason2/CMakeFiles/mason_materializer.dir/flags.make

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o: apps/mason2/CMakeFiles/mason_materializer.dir/flags.make
apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o: /home/svilsen/seqan/seqan-src/apps/mason2/mason_materializer.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/mason2 && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o -c /home/svilsen/seqan/seqan-src/apps/mason2/mason_materializer.cpp

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mason_materializer.dir/mason_materializer.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/mason2 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/mason2/mason_materializer.cpp > CMakeFiles/mason_materializer.dir/mason_materializer.cpp.i

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mason_materializer.dir/mason_materializer.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/mason2 && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/mason2/mason_materializer.cpp -o CMakeFiles/mason_materializer.dir/mason_materializer.cpp.s

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.requires:
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.requires

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.provides: apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.requires
	$(MAKE) -f apps/mason2/CMakeFiles/mason_materializer.dir/build.make apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.provides.build
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.provides

apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.provides.build: apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o

# Object files for target mason_materializer
mason_materializer_OBJECTS = \
"CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o"

# External object files for target mason_materializer
mason_materializer_EXTERNAL_OBJECTS =

bin/mason_materializer: apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o
bin/mason_materializer: apps/mason2/CMakeFiles/mason_materializer.dir/build.make
bin/mason_materializer: /usr/lib/x86_64-linux-gnu/libz.so
bin/mason_materializer: apps/mason2/libmason_sim.a
bin/mason_materializer: /usr/lib/x86_64-linux-gnu/libz.so
bin/mason_materializer: apps/mason2/CMakeFiles/mason_materializer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/mason_materializer"
	cd /home/svilsen/seqan/seqan-build/debug/apps/mason2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mason_materializer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/mason2/CMakeFiles/mason_materializer.dir/build: bin/mason_materializer
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/build

apps/mason2/CMakeFiles/mason_materializer.dir/requires: apps/mason2/CMakeFiles/mason_materializer.dir/mason_materializer.cpp.o.requires
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/requires

apps/mason2/CMakeFiles/mason_materializer.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/mason2 && $(CMAKE_COMMAND) -P CMakeFiles/mason_materializer.dir/cmake_clean.cmake
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/clean

apps/mason2/CMakeFiles/mason_materializer.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/mason2 /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/mason2 /home/svilsen/seqan/seqan-build/debug/apps/mason2/CMakeFiles/mason_materializer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/mason2/CMakeFiles/mason_materializer.dir/depend

