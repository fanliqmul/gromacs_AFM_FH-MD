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
include src/gromacs/fft/tests/CMakeFiles/fft-test.dir/depend.make

# Include the progress variables for this target.
include src/gromacs/fft/tests/CMakeFiles/fft-test.dir/progress.make

# Include the compile flags for this target's objects.
include src/gromacs/fft/tests/CMakeFiles/fft-test.dir/flags.make

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.o: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/flags.make
src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.o: src/gromacs/fft/tests/fft.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fft-test.dir/fft.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests/fft.cpp

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fft-test.dir/fft.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests/fft.cpp > CMakeFiles/fft-test.dir/fft.cpp.i

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fft-test.dir/fft.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests/fft.cpp -o CMakeFiles/fft-test.dir/fft.cpp.s

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/flags.make
src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/testutils/unittest_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp > CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.i

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/testutils/unittest_main.cpp -o CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.s

# Object files for target fft-test
fft__test_OBJECTS = \
"CMakeFiles/fft-test.dir/fft.cpp.o" \
"CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o"

# External object files for target fft-test
fft__test_EXTERNAL_OBJECTS =

bin/fft-test: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/fft.cpp.o
bin/fft-test: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/__/__/__/testutils/unittest_main.cpp.o
bin/fft-test: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/build.make
bin/fft-test: lib/libtestutils.a
bin/fft-test: lib/libgromacs_mpi.so.2.0.0
bin/fft-test: lib/libgmock.a
bin/fft-test: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
bin/fft-test: src/contrib/fftw/fftwBuild-prefix/lib/libfftw3f.a
bin/fft-test: src/gromacs/fft/tests/CMakeFiles/fft-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../../../../bin/fft-test"
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fft-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/gromacs/fft/tests/CMakeFiles/fft-test.dir/build: bin/fft-test

.PHONY : src/gromacs/fft/tests/CMakeFiles/fft-test.dir/build

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/clean:
	cd /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests && $(CMAKE_COMMAND) -P CMakeFiles/fft-test.dir/cmake_clean.cmake
.PHONY : src/gromacs/fft/tests/CMakeFiles/fft-test.dir/clean

src/gromacs/fft/tests/CMakeFiles/fft-test.dir/depend:
	cd /home/fanli/workspace/gromacs_fh_debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests /home/fanli/workspace/gromacs_fh_debug/src/gromacs/fft/tests/CMakeFiles/fft-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/gromacs/fft/tests/CMakeFiles/fft-test.dir/depend

