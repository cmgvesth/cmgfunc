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
CMAKE_SOURCE_DIR = /home/cmgfunc/ProtFun/ProtFun/psipred-2.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cmgfunc/ProtFun/ProtFun/psipred-2.5

# Include any dependencies generated for this target.
include CMakeFiles/./bin/psipass2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/./bin/psipass2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/./bin/psipass2.dir/flags.make

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o: CMakeFiles/./bin/psipass2.dir/flags.make
CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o: src/sspred_hmulti.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cmgfunc/ProtFun/ProtFun/psipred-2.5/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o   -c /home/cmgfunc/ProtFun/ProtFun/psipred-2.5/src/sspred_hmulti.c

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/cmgfunc/ProtFun/ProtFun/psipred-2.5/src/sspred_hmulti.c > CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.i

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/cmgfunc/ProtFun/ProtFun/psipred-2.5/src/sspred_hmulti.c -o CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.s

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.requires:
.PHONY : CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.requires

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.provides: CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.requires
	$(MAKE) -f CMakeFiles/./bin/psipass2.dir/build.make CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.provides.build
.PHONY : CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.provides

CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.provides.build: CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o

# Object files for target ./bin/psipass2
_/bin/psipass2_OBJECTS = \
"CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o"

# External object files for target ./bin/psipass2
_/bin/psipass2_EXTERNAL_OBJECTS =

./bin/psipass2: CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o
./bin/psipass2: CMakeFiles/./bin/psipass2.dir/build.make
./bin/psipass2: CMakeFiles/./bin/psipass2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ./bin/psipass2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/./bin/psipass2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/./bin/psipass2.dir/build: ./bin/psipass2
.PHONY : CMakeFiles/./bin/psipass2.dir/build

CMakeFiles/./bin/psipass2.dir/requires: CMakeFiles/./bin/psipass2.dir/src/sspred_hmulti.c.o.requires
.PHONY : CMakeFiles/./bin/psipass2.dir/requires

CMakeFiles/./bin/psipass2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/./bin/psipass2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/./bin/psipass2.dir/clean

CMakeFiles/./bin/psipass2.dir/depend:
	cd /home/cmgfunc/ProtFun/ProtFun/psipred-2.5 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cmgfunc/ProtFun/ProtFun/psipred-2.5 /home/cmgfunc/ProtFun/ProtFun/psipred-2.5 /home/cmgfunc/ProtFun/ProtFun/psipred-2.5 /home/cmgfunc/ProtFun/ProtFun/psipred-2.5 /home/cmgfunc/ProtFun/ProtFun/psipred-2.5/CMakeFiles/bin/psipass2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/./bin/psipass2.dir/depend
