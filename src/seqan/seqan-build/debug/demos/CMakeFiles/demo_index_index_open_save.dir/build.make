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
include demos/CMakeFiles/demo_index_index_open_save.dir/depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/demo_index_index_open_save.dir/progress.make

# Include the compile flags for this target's objects.
include demos/CMakeFiles/demo_index_index_open_save.dir/flags.make

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o: demos/CMakeFiles/demo_index_index_open_save.dir/flags.make
demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o: /home/svilsen/seqan/seqan-src/demos/index/index_open_save.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o -c /home/svilsen/seqan/seqan-src/demos/index/index_open_save.cpp

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/demos/index/index_open_save.cpp > CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.i

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/demos/index/index_open_save.cpp -o CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.s

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.requires:
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.requires

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.provides: demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.requires
	$(MAKE) -f demos/CMakeFiles/demo_index_index_open_save.dir/build.make demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.provides.build
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.provides

demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.provides.build: demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o

# Object files for target demo_index_index_open_save
demo_index_index_open_save_OBJECTS = \
"CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o"

# External object files for target demo_index_index_open_save
demo_index_index_open_save_EXTERNAL_OBJECTS =

bin/demo_index_index_open_save: demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o
bin/demo_index_index_open_save: demos/CMakeFiles/demo_index_index_open_save.dir/build.make
bin/demo_index_index_open_save: /usr/lib/x86_64-linux-gnu/libz.so
bin/demo_index_index_open_save: /usr/lib/x86_64-linux-gnu/libbz2.so
bin/demo_index_index_open_save: demos/CMakeFiles/demo_index_index_open_save.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/demo_index_index_open_save"
	cd /home/svilsen/seqan/seqan-build/debug/demos && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo_index_index_open_save.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/CMakeFiles/demo_index_index_open_save.dir/build: bin/demo_index_index_open_save
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/build

demos/CMakeFiles/demo_index_index_open_save.dir/requires: demos/CMakeFiles/demo_index_index_open_save.dir/index/index_open_save.cpp.o.requires
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/requires

demos/CMakeFiles/demo_index_index_open_save.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/demos && $(CMAKE_COMMAND) -P CMakeFiles/demo_index_index_open_save.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/clean

demos/CMakeFiles/demo_index_index_open_save.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/demos /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/demos /home/svilsen/seqan/seqan-build/debug/demos/CMakeFiles/demo_index_index_open_save.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/demo_index_index_open_save.dir/depend

