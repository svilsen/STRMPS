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
include apps/samcat/CMakeFiles/samcat.dir/depend.make

# Include the progress variables for this target.
include apps/samcat/CMakeFiles/samcat.dir/progress.make

# Include the compile flags for this target's objects.
include apps/samcat/CMakeFiles/samcat.dir/flags.make

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o: apps/samcat/CMakeFiles/samcat.dir/flags.make
apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o: /home/svilsen/seqan/seqan-src/apps/samcat/samcat.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/samcat && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/samcat.dir/samcat.cpp.o -c /home/svilsen/seqan/seqan-src/apps/samcat/samcat.cpp

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/samcat.dir/samcat.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/samcat && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/samcat/samcat.cpp > CMakeFiles/samcat.dir/samcat.cpp.i

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/samcat.dir/samcat.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/samcat && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/samcat/samcat.cpp -o CMakeFiles/samcat.dir/samcat.cpp.s

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.requires:
.PHONY : apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.requires

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.provides: apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.requires
	$(MAKE) -f apps/samcat/CMakeFiles/samcat.dir/build.make apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.provides.build
.PHONY : apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.provides

apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.provides.build: apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o

# Object files for target samcat
samcat_OBJECTS = \
"CMakeFiles/samcat.dir/samcat.cpp.o"

# External object files for target samcat
samcat_EXTERNAL_OBJECTS =

bin/samcat: apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o
bin/samcat: apps/samcat/CMakeFiles/samcat.dir/build.make
bin/samcat: /usr/lib/x86_64-linux-gnu/libz.so
bin/samcat: apps/samcat/CMakeFiles/samcat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/samcat"
	cd /home/svilsen/seqan/seqan-build/debug/apps/samcat && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/samcat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/samcat/CMakeFiles/samcat.dir/build: bin/samcat
.PHONY : apps/samcat/CMakeFiles/samcat.dir/build

apps/samcat/CMakeFiles/samcat.dir/requires: apps/samcat/CMakeFiles/samcat.dir/samcat.cpp.o.requires
.PHONY : apps/samcat/CMakeFiles/samcat.dir/requires

apps/samcat/CMakeFiles/samcat.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/samcat && $(CMAKE_COMMAND) -P CMakeFiles/samcat.dir/cmake_clean.cmake
.PHONY : apps/samcat/CMakeFiles/samcat.dir/clean

apps/samcat/CMakeFiles/samcat.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/samcat /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/samcat /home/svilsen/seqan/seqan-build/debug/apps/samcat/CMakeFiles/samcat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/samcat/CMakeFiles/samcat.dir/depend
