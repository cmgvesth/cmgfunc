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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin

# Include any dependencies generated for this target.
include CMakeFiles/pseg.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/pseg.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pseg.dir/flags.make

CMakeFiles/pseg.dir/pseg/pseg.c.o: CMakeFiles/pseg.dir/flags.make
CMakeFiles/pseg.dir/pseg/pseg.c.o: ../pseg/pseg.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/pseg.dir/pseg/pseg.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/pseg.dir/pseg/pseg.c.o   -c /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pseg.c

CMakeFiles/pseg.dir/pseg/pseg.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pseg.dir/pseg/pseg.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pseg.c > CMakeFiles/pseg.dir/pseg/pseg.c.i

CMakeFiles/pseg.dir/pseg/pseg.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pseg.dir/pseg/pseg.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pseg.c -o CMakeFiles/pseg.dir/pseg/pseg.c.s

CMakeFiles/pseg.dir/pseg/pseg.c.o.requires:
.PHONY : CMakeFiles/pseg.dir/pseg/pseg.c.o.requires

CMakeFiles/pseg.dir/pseg/pseg.c.o.provides: CMakeFiles/pseg.dir/pseg/pseg.c.o.requires
	$(MAKE) -f CMakeFiles/pseg.dir/build.make CMakeFiles/pseg.dir/pseg/pseg.c.o.provides.build
.PHONY : CMakeFiles/pseg.dir/pseg/pseg.c.o.provides

CMakeFiles/pseg.dir/pseg/pseg.c.o.provides.build: CMakeFiles/pseg.dir/pseg/pseg.c.o

CMakeFiles/pseg.dir/pseg/pgenwin.c.o: CMakeFiles/pseg.dir/flags.make
CMakeFiles/pseg.dir/pseg/pgenwin.c.o: ../pseg/pgenwin.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/pseg.dir/pseg/pgenwin.c.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/pseg.dir/pseg/pgenwin.c.o   -c /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pgenwin.c

CMakeFiles/pseg.dir/pseg/pgenwin.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/pseg.dir/pseg/pgenwin.c.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pgenwin.c > CMakeFiles/pseg.dir/pseg/pgenwin.c.i

CMakeFiles/pseg.dir/pseg/pgenwin.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/pseg.dir/pseg/pgenwin.c.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/pseg/pgenwin.c -o CMakeFiles/pseg.dir/pseg/pgenwin.c.s

CMakeFiles/pseg.dir/pseg/pgenwin.c.o.requires:
.PHONY : CMakeFiles/pseg.dir/pseg/pgenwin.c.o.requires

CMakeFiles/pseg.dir/pseg/pgenwin.c.o.provides: CMakeFiles/pseg.dir/pseg/pgenwin.c.o.requires
	$(MAKE) -f CMakeFiles/pseg.dir/build.make CMakeFiles/pseg.dir/pseg/pgenwin.c.o.provides.build
.PHONY : CMakeFiles/pseg.dir/pseg/pgenwin.c.o.provides

CMakeFiles/pseg.dir/pseg/pgenwin.c.o.provides.build: CMakeFiles/pseg.dir/pseg/pgenwin.c.o

# Object files for target pseg
pseg_OBJECTS = \
"CMakeFiles/pseg.dir/pseg/pseg.c.o" \
"CMakeFiles/pseg.dir/pseg/pgenwin.c.o"

# External object files for target pseg
pseg_EXTERNAL_OBJECTS =

pseg: CMakeFiles/pseg.dir/pseg/pseg.c.o
pseg: CMakeFiles/pseg.dir/pseg/pgenwin.c.o
pseg: CMakeFiles/pseg.dir/build.make
pseg: CMakeFiles/pseg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable pseg"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pseg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pseg.dir/build: pseg
.PHONY : CMakeFiles/pseg.dir/build

CMakeFiles/pseg.dir/requires: CMakeFiles/pseg.dir/pseg/pseg.c.o.requires
CMakeFiles/pseg.dir/requires: CMakeFiles/pseg.dir/pseg/pgenwin.c.o.requires
.PHONY : CMakeFiles/pseg.dir/requires

CMakeFiles/pseg.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pseg.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pseg.dir/clean

CMakeFiles/pseg.dir/depend:
	cd /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1 /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1 /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin /home/steven/Projects/CMG-Biotools-Share/protfun_cmake/seg-1/bin/CMakeFiles/pseg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pseg.dir/depend

