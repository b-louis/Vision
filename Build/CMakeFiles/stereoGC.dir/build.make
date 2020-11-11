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
CMAKE_SOURCE_DIR = /home/x1/M2/VISION/gcdispar

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/x1/M2/VISION/gcdispar/Build

# Include any dependencies generated for this target.
include CMakeFiles/stereoGC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/stereoGC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/stereoGC.dir/flags.make

CMakeFiles/stereoGC.dir/stereoGC.cpp.o: CMakeFiles/stereoGC.dir/flags.make
CMakeFiles/stereoGC.dir/stereoGC.cpp.o: ../stereoGC.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/x1/M2/VISION/gcdispar/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/stereoGC.dir/stereoGC.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stereoGC.dir/stereoGC.cpp.o -c /home/x1/M2/VISION/gcdispar/stereoGC.cpp

CMakeFiles/stereoGC.dir/stereoGC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stereoGC.dir/stereoGC.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/x1/M2/VISION/gcdispar/stereoGC.cpp > CMakeFiles/stereoGC.dir/stereoGC.cpp.i

CMakeFiles/stereoGC.dir/stereoGC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stereoGC.dir/stereoGC.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/x1/M2/VISION/gcdispar/stereoGC.cpp -o CMakeFiles/stereoGC.dir/stereoGC.cpp.s

CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o: CMakeFiles/stereoGC.dir/flags.make
CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o: ../maxflow/graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/x1/M2/VISION/gcdispar/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o -c /home/x1/M2/VISION/gcdispar/maxflow/graph.cpp

CMakeFiles/stereoGC.dir/maxflow/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/stereoGC.dir/maxflow/graph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/x1/M2/VISION/gcdispar/maxflow/graph.cpp > CMakeFiles/stereoGC.dir/maxflow/graph.cpp.i

CMakeFiles/stereoGC.dir/maxflow/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/stereoGC.dir/maxflow/graph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/x1/M2/VISION/gcdispar/maxflow/graph.cpp -o CMakeFiles/stereoGC.dir/maxflow/graph.cpp.s

# Object files for target stereoGC
stereoGC_OBJECTS = \
"CMakeFiles/stereoGC.dir/stereoGC.cpp.o" \
"CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o"

# External object files for target stereoGC
stereoGC_EXTERNAL_OBJECTS =

stereoGC: CMakeFiles/stereoGC.dir/stereoGC.cpp.o
stereoGC: CMakeFiles/stereoGC.dir/maxflow/graph.cpp.o
stereoGC: CMakeFiles/stereoGC.dir/build.make
stereoGC: /home/x1/anaconda3/lib/libQt5OpenGL.so.5.9.7
stereoGC: /usr/lib/x86_64-linux-gnu/libGL.so
stereoGC: /usr/lib/x86_64-linux-gnu/libGLU.so
stereoGC: /home/x1/anaconda3/lib/libQt5Widgets.so.5.9.7
stereoGC: /home/x1/anaconda3/lib/libQt5Gui.so.5.9.7
stereoGC: /home/x1/anaconda3/lib/libQt5Core.so.5.9.7
stereoGC: CMakeFiles/stereoGC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/x1/M2/VISION/gcdispar/Build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable stereoGC"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/stereoGC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/stereoGC.dir/build: stereoGC

.PHONY : CMakeFiles/stereoGC.dir/build

CMakeFiles/stereoGC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/stereoGC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/stereoGC.dir/clean

CMakeFiles/stereoGC.dir/depend:
	cd /home/x1/M2/VISION/gcdispar/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/x1/M2/VISION/gcdispar /home/x1/M2/VISION/gcdispar /home/x1/M2/VISION/gcdispar/Build /home/x1/M2/VISION/gcdispar/Build /home/x1/M2/VISION/gcdispar/Build/CMakeFiles/stereoGC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/stereoGC.dir/depend

