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
include demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/depend.make

# Include the progress variables for this target.
include demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/progress.make

# Include the compile flags for this target's objects.
include demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/flags.make

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/flags.make
demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o: /home/svilsen/seqan/seqan-src/demos/tutorial/rnaseq/genequant_assignment3.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o -c /home/svilsen/seqan/seqan-src/demos/tutorial/rnaseq/genequant_assignment3.cpp

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/demos/tutorial/rnaseq/genequant_assignment3.cpp > CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.i

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/demos && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/demos/tutorial/rnaseq/genequant_assignment3.cpp -o CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.s

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.requires:
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.requires

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.provides: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.requires
	$(MAKE) -f demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/build.make demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.provides.build
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.provides

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.provides.build: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o

# Object files for target demo_tutorial_rnaseq_genequant_assignment3
demo_tutorial_rnaseq_genequant_assignment3_OBJECTS = \
"CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o"

# External object files for target demo_tutorial_rnaseq_genequant_assignment3
demo_tutorial_rnaseq_genequant_assignment3_EXTERNAL_OBJECTS =

bin/demo_tutorial_rnaseq_genequant_assignment3: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o
bin/demo_tutorial_rnaseq_genequant_assignment3: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/build.make
bin/demo_tutorial_rnaseq_genequant_assignment3: /usr/lib/x86_64-linux-gnu/libz.so
bin/demo_tutorial_rnaseq_genequant_assignment3: /usr/lib/x86_64-linux-gnu/libbz2.so
bin/demo_tutorial_rnaseq_genequant_assignment3: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/demo_tutorial_rnaseq_genequant_assignment3"
	cd /home/svilsen/seqan/seqan-build/debug/demos && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/build: bin/demo_tutorial_rnaseq_genequant_assignment3
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/build

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/requires: demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/tutorial/rnaseq/genequant_assignment3.cpp.o.requires
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/requires

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/demos && $(CMAKE_COMMAND) -P CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/cmake_clean.cmake
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/clean

demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/demos /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/demos /home/svilsen/seqan/seqan-build/debug/demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : demos/CMakeFiles/demo_tutorial_rnaseq_genequant_assignment3.dir/depend

