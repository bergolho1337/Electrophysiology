# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_SOURCE_DIR = /home/berg/Github/Electrophysiology/Fenton-Experiment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/berg/Github/Electrophysiology/Fenton-Experiment/build

# Include any dependencies generated for this target.
include src/config/CMakeFiles/config.dir/depend.make

# Include the progress variables for this target.
include src/config/CMakeFiles/config.dir/progress.make

# Include the compile flags for this target's objects.
include src/config/CMakeFiles/config.dir/flags.make

src/config/CMakeFiles/config.dir/config.cpp.o: src/config/CMakeFiles/config.dir/flags.make
src/config/CMakeFiles/config.dir/config.cpp.o: ../src/config/config.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/berg/Github/Electrophysiology/Fenton-Experiment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/config/CMakeFiles/config.dir/config.cpp.o"
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/config.dir/config.cpp.o -c /home/berg/Github/Electrophysiology/Fenton-Experiment/src/config/config.cpp

src/config/CMakeFiles/config.dir/config.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/config.dir/config.cpp.i"
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/berg/Github/Electrophysiology/Fenton-Experiment/src/config/config.cpp > CMakeFiles/config.dir/config.cpp.i

src/config/CMakeFiles/config.dir/config.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/config.dir/config.cpp.s"
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/berg/Github/Electrophysiology/Fenton-Experiment/src/config/config.cpp -o CMakeFiles/config.dir/config.cpp.s

# Object files for target config
config_OBJECTS = \
"CMakeFiles/config.dir/config.cpp.o"

# External object files for target config
config_EXTERNAL_OBJECTS =

src/config/libconfig.a: src/config/CMakeFiles/config.dir/config.cpp.o
src/config/libconfig.a: src/config/CMakeFiles/config.dir/build.make
src/config/libconfig.a: src/config/CMakeFiles/config.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/berg/Github/Electrophysiology/Fenton-Experiment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libconfig.a"
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && $(CMAKE_COMMAND) -P CMakeFiles/config.dir/cmake_clean_target.cmake
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/config.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/config/CMakeFiles/config.dir/build: src/config/libconfig.a

.PHONY : src/config/CMakeFiles/config.dir/build

src/config/CMakeFiles/config.dir/clean:
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config && $(CMAKE_COMMAND) -P CMakeFiles/config.dir/cmake_clean.cmake
.PHONY : src/config/CMakeFiles/config.dir/clean

src/config/CMakeFiles/config.dir/depend:
	cd /home/berg/Github/Electrophysiology/Fenton-Experiment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/berg/Github/Electrophysiology/Fenton-Experiment /home/berg/Github/Electrophysiology/Fenton-Experiment/src/config /home/berg/Github/Electrophysiology/Fenton-Experiment/build /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config /home/berg/Github/Electrophysiology/Fenton-Experiment/build/src/config/CMakeFiles/config.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/config/CMakeFiles/config.dir/depend

