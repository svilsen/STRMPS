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
include apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/depend.make

# Include the progress variables for this target.
include apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/progress.make

# Include the compile flags for this target's objects.
include apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/flags.make

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/flags.make
apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o: /home/svilsen/seqan/seqan-src/apps/param_chooser/test_param_chooser.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/svilsen/seqan/seqan-build/debug/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o"
	cd /home/svilsen/seqan/seqan-build/debug/apps/param_chooser && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o -c /home/svilsen/seqan/seqan-src/apps/param_chooser/test_param_chooser.cpp

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.i"
	cd /home/svilsen/seqan/seqan-build/debug/apps/param_chooser && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/svilsen/seqan/seqan-src/apps/param_chooser/test_param_chooser.cpp > CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.i

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.s"
	cd /home/svilsen/seqan/seqan-build/debug/apps/param_chooser && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/svilsen/seqan/seqan-src/apps/param_chooser/test_param_chooser.cpp -o CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.s

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.requires:
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.requires

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.provides: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.requires
	$(MAKE) -f apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/build.make apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.provides.build
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.provides

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.provides.build: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o

# Object files for target test_funcs_param_chooser
test_funcs_param_chooser_OBJECTS = \
"CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o"

# External object files for target test_funcs_param_chooser
test_funcs_param_chooser_EXTERNAL_OBJECTS =

bin/test_funcs_param_chooser: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o
bin/test_funcs_param_chooser: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/build.make
bin/test_funcs_param_chooser: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../bin/test_funcs_param_chooser"
	cd /home/svilsen/seqan/seqan-build/debug/apps/param_chooser && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_funcs_param_chooser.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/build: bin/test_funcs_param_chooser
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/build

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/requires: apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/test_param_chooser.cpp.o.requires
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/requires

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/clean:
	cd /home/svilsen/seqan/seqan-build/debug/apps/param_chooser && $(CMAKE_COMMAND) -P CMakeFiles/test_funcs_param_chooser.dir/cmake_clean.cmake
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/clean

apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/depend:
	cd /home/svilsen/seqan/seqan-build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/svilsen/seqan/seqan-src /home/svilsen/seqan/seqan-src/apps/param_chooser /home/svilsen/seqan/seqan-build/debug /home/svilsen/seqan/seqan-build/debug/apps/param_chooser /home/svilsen/seqan/seqan-build/debug/apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/param_chooser/CMakeFiles/test_funcs_param_chooser.dir/depend
