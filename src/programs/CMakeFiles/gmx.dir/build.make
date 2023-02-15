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
include src/programs/CMakeFiles/gmx.dir/depend.make

# Include the progress variables for this target.
include src/programs/CMakeFiles/gmx.dir/progress.make

# Include the compile flags for this target's objects.
include src/programs/CMakeFiles/gmx.dir/flags.make

src/programs/CMakeFiles/gmx.dir/gmx.cpp.o: src/programs/CMakeFiles/gmx.dir/flags.make
src/programs/CMakeFiles/gmx.dir/gmx.cpp.o: src/programs/gmx.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/programs/CMakeFiles/gmx.dir/gmx.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmx.dir/gmx.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/programs/gmx.cpp

src/programs/CMakeFiles/gmx.dir/gmx.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmx.dir/gmx.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/programs/gmx.cpp > CMakeFiles/gmx.dir/gmx.cpp.i

src/programs/CMakeFiles/gmx.dir/gmx.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmx.dir/gmx.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/programs/gmx.cpp -o CMakeFiles/gmx.dir/gmx.cpp.s

src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.o: src/programs/CMakeFiles/gmx.dir/flags.make
src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.o: src/programs/legacymodules.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.o"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmx.dir/legacymodules.cpp.o -c /home/fanli/workspace/gromacs_fh_debug/src/programs/legacymodules.cpp

src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmx.dir/legacymodules.cpp.i"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/fanli/workspace/gromacs_fh_debug/src/programs/legacymodules.cpp > CMakeFiles/gmx.dir/legacymodules.cpp.i

src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmx.dir/legacymodules.cpp.s"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/fanli/workspace/gromacs_fh_debug/src/programs/legacymodules.cpp -o CMakeFiles/gmx.dir/legacymodules.cpp.s

# Object files for target gmx
gmx_OBJECTS = \
"CMakeFiles/gmx.dir/gmx.cpp.o" \
"CMakeFiles/gmx.dir/legacymodules.cpp.o"

# External object files for target gmx
gmx_EXTERNAL_OBJECTS = \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/coupling.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/estimate.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/fh.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/flow_velocity.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/init.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/interpolation.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_md_support.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_stat.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_tgroup.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_update.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/output.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/parser.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/sfunction.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/md.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/mdrun.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/membed.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/repl_ex.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/resource-division.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/runner.cpp.o" \
"/home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/view_objlib.dir/view/view.cpp.o"

bin/gmx_mpi: src/programs/CMakeFiles/gmx.dir/gmx.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/gmx.dir/legacymodules.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/coupling.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/estimate.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/fh.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/flow_velocity.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/init.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/interpolation.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_md_support.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_stat.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_tgroup.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/new_update.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/output.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/parser.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/__/gromacs/fhmdlib/sfunction.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/md.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/mdrun.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/membed.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/repl_ex.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/resource-division.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/mdrun_objlib.dir/mdrun/runner.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/view_objlib.dir/view/view.cpp.o
bin/gmx_mpi: src/programs/CMakeFiles/gmx.dir/build.make
bin/gmx_mpi: lib/libgromacs_mpi.so.2.0.0
bin/gmx_mpi: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
bin/gmx_mpi: src/contrib/fftw/fftwBuild-prefix/lib/libfftw3f.a
bin/gmx_mpi: src/programs/CMakeFiles/gmx.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/fanli/workspace/gromacs_fh_debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../../bin/gmx_mpi"
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmx.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/programs/CMakeFiles/gmx.dir/build: bin/gmx_mpi

.PHONY : src/programs/CMakeFiles/gmx.dir/build

src/programs/CMakeFiles/gmx.dir/clean:
	cd /home/fanli/workspace/gromacs_fh_debug/src/programs && $(CMAKE_COMMAND) -P CMakeFiles/gmx.dir/cmake_clean.cmake
.PHONY : src/programs/CMakeFiles/gmx.dir/clean

src/programs/CMakeFiles/gmx.dir/depend:
	cd /home/fanli/workspace/gromacs_fh_debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/programs /home/fanli/workspace/gromacs_fh_debug /home/fanli/workspace/gromacs_fh_debug/src/programs /home/fanli/workspace/gromacs_fh_debug/src/programs/CMakeFiles/gmx.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/programs/CMakeFiles/gmx.dir/depend
