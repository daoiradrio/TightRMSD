# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/baum/TightRMSD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/baum/TightRMSD/build

# Include any dependencies generated for this target.
include CMakeFiles/rmsd.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/rmsd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rmsd.dir/flags.make

CMakeFiles/rmsd.dir/main_files/main.cpp.o: CMakeFiles/rmsd.dir/flags.make
CMakeFiles/rmsd.dir/main_files/main.cpp.o: ../main_files/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/baum/TightRMSD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rmsd.dir/main_files/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rmsd.dir/main_files/main.cpp.o -c /home/baum/TightRMSD/main_files/main.cpp

CMakeFiles/rmsd.dir/main_files/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rmsd.dir/main_files/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/baum/TightRMSD/main_files/main.cpp > CMakeFiles/rmsd.dir/main_files/main.cpp.i

CMakeFiles/rmsd.dir/main_files/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rmsd.dir/main_files/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/baum/TightRMSD/main_files/main.cpp -o CMakeFiles/rmsd.dir/main_files/main.cpp.s

# Object files for target rmsd
rmsd_OBJECTS = \
"CMakeFiles/rmsd.dir/main_files/main.cpp.o"

# External object files for target rmsd
rmsd_EXTERNAL_OBJECTS =

rmsd: CMakeFiles/rmsd.dir/main_files/main.cpp.o
rmsd: CMakeFiles/rmsd.dir/build.make
rmsd: CMakeFiles/rmsd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/baum/TightRMSD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable rmsd"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rmsd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/rmsd.dir/build: rmsd

.PHONY : CMakeFiles/rmsd.dir/build

CMakeFiles/rmsd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rmsd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rmsd.dir/clean

CMakeFiles/rmsd.dir/depend:
	cd /home/baum/TightRMSD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/baum/TightRMSD /home/baum/TightRMSD /home/baum/TightRMSD/build /home/baum/TightRMSD/build /home/baum/TightRMSD/build/CMakeFiles/rmsd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rmsd.dir/depend

