# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/fanli/workspace/gromacs_fh_debug

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fanli/workspace/gromacs_fh_debug

# Include any dependencies generated for this target.
include src/gromacs/simd/tests/CMakeFiles/simd-test.dir/depend.make

# Include the progress variables for this target.
include src/gromacs/simd/tests/CMakeFiles/simd-test.dir/progress.make

# Include the compile flags for this target's objects.
include src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o: src/gromacs/simd/tests/bootstrap_loadstore.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/bootstrap_loadstore.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/bootstrap_loadstore.cpp > CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/bootstrap_loadstore.cpp -o CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.o: src/gromacs/simd/tests/base.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/base.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/base.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/base.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/base.cpp > CMakeFiles/simd-test.dir/base.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/base.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/base.cpp -o CMakeFiles/simd-test.dir/base.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.o: src/gromacs/simd/tests/simd.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd.cpp > CMakeFiles/simd-test.dir/simd.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd.cpp -o CMakeFiles/simd-test.dir/simd.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o: src/gromacs/simd/tests/simd_floatingpoint.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint.cpp > CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint.cpp -o CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o: src/gromacs/simd/tests/simd_floatingpoint_util.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint_util.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint_util.cpp > CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_floatingpoint_util.cpp -o CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o: src/gromacs/simd/tests/simd_vector_operations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_vector_operations.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd_vector_operations.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_vector_operations.cpp > CMakeFiles/simd-test.dir/simd_vector_operations.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd_vector_operations.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_vector_operations.cpp -o CMakeFiles/simd-test.dir/simd_vector_operations.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.o: src/gromacs/simd/tests/simd_math.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd_math.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_math.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd_math.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_math.cpp > CMakeFiles/simd-test.dir/simd_math.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd_math.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_math.cpp -o CMakeFiles/simd-test.dir/simd_math.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.o: src/gromacs/simd/tests/simd_integer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd_integer.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_integer.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd_integer.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_integer.cpp > CMakeFiles/simd-test.dir/simd_integer.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd_integer.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd_integer.cpp -o CMakeFiles/simd-test.dir/simd_integer.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.o: src/gromacs/simd/tests/simd4.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd4.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd4.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4.cpp > CMakeFiles/simd-test.dir/simd4.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd4.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4.cpp -o CMakeFiles/simd-test.dir/simd4.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o: src/gromacs/simd/tests/simd4_floatingpoint.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_floatingpoint.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_floatingpoint.cpp > CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_floatingpoint.cpp -o CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o: src/gromacs/simd/tests/simd4_vector_operations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_vector_operations.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_vector_operations.cpp > CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_vector_operations.cpp -o CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.o: src/gromacs/simd/tests/simd4_math.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/simd4_math.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_math.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/simd4_math.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_math.cpp > CMakeFiles/simd-test.dir/simd4_math.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/simd4_math.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/simd4_math.cpp -o CMakeFiles/simd-test.dir/simd4_math.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.o: src/gromacs/simd/tests/scalar.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/scalar.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/scalar.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar.cpp > CMakeFiles/simd-test.dir/scalar.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/scalar.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar.cpp -o CMakeFiles/simd-test.dir/scalar.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.o: src/gromacs/simd/tests/scalar_util.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/scalar_util.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_util.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/scalar_util.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_util.cpp > CMakeFiles/simd-test.dir/scalar_util.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/scalar_util.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_util.cpp -o CMakeFiles/simd-test.dir/scalar_util.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.o: src/gromacs/simd/tests/scalar_math.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/scalar_math.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_math.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/scalar_math.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_math.cpp > CMakeFiles/simd-test.dir/scalar_math.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/scalar_math.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/scalar_math.cpp -o CMakeFiles/simd-test.dir/scalar_math.cpp.s

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/flags.make
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/testutils/unittest_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp > CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.i

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp -o CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.s

# Object files for target simd-test
simd__test_OBJECTS = \
"CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o" \
"CMakeFiles/simd-test.dir/base.cpp.o" \
"CMakeFiles/simd-test.dir/simd.cpp.o" \
"CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o" \
"CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o" \
"CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o" \
"CMakeFiles/simd-test.dir/simd_math.cpp.o" \
"CMakeFiles/simd-test.dir/simd_integer.cpp.o" \
"CMakeFiles/simd-test.dir/simd4.cpp.o" \
"CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o" \
"CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o" \
"CMakeFiles/simd-test.dir/simd4_math.cpp.o" \
"CMakeFiles/simd-test.dir/scalar.cpp.o" \
"CMakeFiles/simd-test.dir/scalar_util.cpp.o" \
"CMakeFiles/simd-test.dir/scalar_math.cpp.o" \
"CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o"

# External object files for target simd-test
simd__test_EXTERNAL_OBJECTS =

bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/bootstrap_loadstore.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/base.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_floatingpoint_util.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_vector_operations.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_math.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd_integer.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_floatingpoint.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_vector_operations.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/simd4_math.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_util.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/scalar_math.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/__/__/__/testutils/unittest_main.cpp.o
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/build.make
bin/simd-test: lib/libtestutils.a
bin/simd-test: lib/libgromacs_mpi.so.2.0.0
bin/simd-test: lib/libgmock.a
bin/simd-test: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
bin/simd-test: src/contrib/fftw/fftwBuild-prefix/lib/libfftw3f.a
bin/simd-test: src/gromacs/simd/tests/CMakeFiles/simd-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Linking CXX executable ../../../../bin/simd-test"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simd-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/gromacs/simd/tests/CMakeFiles/simd-test.dir/build: bin/simd-test

.PHONY : src/gromacs/simd/tests/CMakeFiles/simd-test.dir/build

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/clean:
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests && $(CMAKE_COMMAND) -P CMakeFiles/simd-test.dir/cmake_clean.cmake
.PHONY : src/gromacs/simd/tests/CMakeFiles/simd-test.dir/clean

src/gromacs/simd/tests/CMakeFiles/simd-test.dir/depend:
	cd /home/fanli/workspace/gromacs_fh_debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests /home/fanli/workspace/gromacs_fh_debug/src/gromacs/simd/tests/CMakeFiles/simd-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/gromacs/simd/tests/CMakeFiles/simd-test.dir/depend

