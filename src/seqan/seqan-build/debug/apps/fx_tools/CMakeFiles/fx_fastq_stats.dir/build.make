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
include apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/depend.make

# Include the progress variables for this target.
include apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/progress.make

# Include the compile flags for this target's objects.
include apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/flags.make

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/flags.make
apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o: /home/svilsen/seqan/seqan-src/apps/fx_tools/fx_fastq_stats.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fx_tools && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o -c /home/svilsen/seqan/seqan-src/apps/fx_tools/fx_fastq_stats.cpp

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fx_tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/fx_tools/fx_fastq_stats.cpp > CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.i

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fx_tools && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/fx_tools/fx_fastq_stats.cpp -o CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.s

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.requires:
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.requires

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.provides: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.requires
	$(MAKE) -f apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/build.make apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.provides.build
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.provides

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.provides.build: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o

# Object files for target fx_fastq_stats
fx_fastq_stats_OBJECTS = \
"CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o"

# External object files for target fx_fastq_stats
fx_fastq_stats_EXTERNAL_OBJECTS =

bin/fx_fastq_stats: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o
bin/fx_fastq_stats: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/build.make
bin/fx_fastq_stats: /usr/lib/x86_64-linux-gnu/libz.so
bin/fx_fastq_stats: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/fx_fastq_stats"
	cd /home/svilsen/seqan/seqan-build/debug/apps/fx_tools && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fx_fastq_stats.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/build: bin/fx_fastq_stats
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/build

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/requires: apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/fx_fastq_stats.cpp.o.requires
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/requires

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/fx_tools && $(CMAKE_COMMAND) -P CMakeFiles/fx_fastq_stats.dir/cmake_clean.cmake
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/clean

apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/fx_tools /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/fx_tools /home/svilsen/seqan/seqan-build/debug/apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/fx_tools/CMakeFiles/fx_fastq_stats.dir/depend

