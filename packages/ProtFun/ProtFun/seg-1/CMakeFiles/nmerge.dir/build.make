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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cmgfunc/ProtFun/ProtFun/seg-1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cmgfunc/ProtFun/ProtFun/seg-1

# Include any dependencies generated for this target.
include CMakeFiles/nmerge.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/nmerge.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nmerge.dir/flags.make

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o: CMakeFiles/nmerge.dir/flags.make
CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o: src/nseg/nmerge.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmgfunc/ProtFun/ProtFun/seg-1/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o   -c /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/nmerge.c

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/nmerge.dir/src/nseg/nmerge.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/nmerge.c > CMakeFiles/nmerge.dir/src/nseg/nmerge.c.i

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/nmerge.dir/src/nseg/nmerge.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/nmerge.c -o CMakeFiles/nmerge.dir/src/nseg/nmerge.c.s

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.requires:
.PHONY : CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.requires

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.provides: CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.requires
	$(MAKE) -f CMakeFiles/nmerge.dir/build.make CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.provides.build
.PHONY : CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.provides

CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.provides.build: CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o

CMakeFiles/nmerge.dir/src/nseg/genwin.c.o: CMakeFiles/nmerge.dir/flags.make
CMakeFiles/nmerge.dir/src/nseg/genwin.c.o: src/nseg/genwin.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmgfunc/ProtFun/ProtFun/seg-1/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/nmerge.dir/src/nseg/genwin.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/nmerge.dir/src/nseg/genwin.c.o   -c /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/genwin.c

CMakeFiles/nmerge.dir/src/nseg/genwin.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/nmerge.dir/src/nseg/genwin.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/genwin.c > CMakeFiles/nmerge.dir/src/nseg/genwin.c.i

CMakeFiles/nmerge.dir/src/nseg/genwin.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/nmerge.dir/src/nseg/genwin.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cmgfunc/ProtFun/ProtFun/seg-1/src/nseg/genwin.c -o CMakeFiles/nmerge.dir/src/nseg/genwin.c.s

CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.requires:
.PHONY : CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.requires

CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.provides: CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.requires
	$(MAKE) -f CMakeFiles/nmerge.dir/build.make CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.provides.build
.PHONY : CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.provides

CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.provides.build: CMakeFiles/nmerge.dir/src/nseg/genwin.c.o

# Object files for target nmerge
nmerge_OBJECTS = \
"CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o" \
"CMakeFiles/nmerge.dir/src/nseg/genwin.c.o"

# External object files for target nmerge
nmerge_EXTERNAL_OBJECTS =

nmerge: CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o
nmerge: CMakeFiles/nmerge.dir/src/nseg/genwin.c.o
nmerge: CMakeFiles/nmerge.dir/build.make
nmerge: CMakeFiles/nmerge.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable nmerge"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nmerge.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/nmerge.dir/build: nmerge
.PHONY : CMakeFiles/nmerge.dir/build

CMakeFiles/nmerge.dir/requires: CMakeFiles/nmerge.dir/src/nseg/nmerge.c.o.requires
CMakeFiles/nmerge.dir/requires: CMakeFiles/nmerge.dir/src/nseg/genwin.c.o.requires
.PHONY : CMakeFiles/nmerge.dir/requires

CMakeFiles/nmerge.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nmerge.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nmerge.dir/clean

CMakeFiles/nmerge.dir/depend:
	cd /home/cmgfunc/ProtFun/ProtFun/seg-1 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmgfunc/ProtFun/ProtFun/seg-1 /home/cmgfunc/ProtFun/ProtFun/seg-1 /home/cmgfunc/ProtFun/ProtFun/seg-1 /home/cmgfunc/ProtFun/ProtFun/seg-1 /home/cmgfunc/ProtFun/ProtFun/seg-1/CMakeFiles/nmerge.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nmerge.dir/depend

