# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /opt/GoogleCode/popot/branches/popot-1.20

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /opt/GoogleCode/popot/branches/popot-1.20/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/example-001.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/example-001.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/example-001.dir/flags.make

examples/CMakeFiles/example-001.dir/example-001.cc.o: examples/CMakeFiles/example-001.dir/flags.make
examples/CMakeFiles/example-001.dir/example-001.cc.o: ../examples/example-001.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /opt/GoogleCode/popot/branches/popot-1.20/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object examples/CMakeFiles/example-001.dir/example-001.cc.o"
	cd /opt/GoogleCode/popot/branches/popot-1.20/build/examples && /usr/lib/ccache/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/example-001.dir/example-001.cc.o -c /opt/GoogleCode/popot/branches/popot-1.20/examples/example-001.cc

examples/CMakeFiles/example-001.dir/example-001.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/example-001.dir/example-001.cc.i"
	cd /opt/GoogleCode/popot/branches/popot-1.20/build/examples && /usr/lib/ccache/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /opt/GoogleCode/popot/branches/popot-1.20/examples/example-001.cc > CMakeFiles/example-001.dir/example-001.cc.i

examples/CMakeFiles/example-001.dir/example-001.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/example-001.dir/example-001.cc.s"
	cd /opt/GoogleCode/popot/branches/popot-1.20/build/examples && /usr/lib/ccache/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /opt/GoogleCode/popot/branches/popot-1.20/examples/example-001.cc -o CMakeFiles/example-001.dir/example-001.cc.s

examples/CMakeFiles/example-001.dir/example-001.cc.o.requires:
.PHONY : examples/CMakeFiles/example-001.dir/example-001.cc.o.requires

examples/CMakeFiles/example-001.dir/example-001.cc.o.provides: examples/CMakeFiles/example-001.dir/example-001.cc.o.requires
	$(MAKE) -f examples/CMakeFiles/example-001.dir/build.make examples/CMakeFiles/example-001.dir/example-001.cc.o.provides.build
.PHONY : examples/CMakeFiles/example-001.dir/example-001.cc.o.provides

examples/CMakeFiles/example-001.dir/example-001.cc.o.provides.build: examples/CMakeFiles/example-001.dir/example-001.cc.o

# Object files for target example-001
example__001_OBJECTS = \
"CMakeFiles/example-001.dir/example-001.cc.o"

# External object files for target example-001
example__001_EXTERNAL_OBJECTS =

examples/example-001: examples/CMakeFiles/example-001.dir/example-001.cc.o
examples/example-001: src/libpopot.so
examples/example-001: examples/CMakeFiles/example-001.dir/build.make
examples/example-001: examples/CMakeFiles/example-001.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable example-001"
	cd /opt/GoogleCode/popot/branches/popot-1.20/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example-001.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/example-001.dir/build: examples/example-001
.PHONY : examples/CMakeFiles/example-001.dir/build

examples/CMakeFiles/example-001.dir/requires: examples/CMakeFiles/example-001.dir/example-001.cc.o.requires
.PHONY : examples/CMakeFiles/example-001.dir/requires

examples/CMakeFiles/example-001.dir/clean:
	cd /opt/GoogleCode/popot/branches/popot-1.20/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/example-001.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/example-001.dir/clean

examples/CMakeFiles/example-001.dir/depend:
	cd /opt/GoogleCode/popot/branches/popot-1.20/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/GoogleCode/popot/branches/popot-1.20 /opt/GoogleCode/popot/branches/popot-1.20/examples /opt/GoogleCode/popot/branches/popot-1.20/build /opt/GoogleCode/popot/branches/popot-1.20/build/examples /opt/GoogleCode/popot/branches/popot-1.20/build/examples/CMakeFiles/example-001.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/example-001.dir/depend

