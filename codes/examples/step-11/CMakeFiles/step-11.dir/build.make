# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.1

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
CMAKE_COMMAND = /opt/apps/cmake/3.1.3/bin/cmake

# The command to remove a file.
RM = /opt/apps/cmake/3.1.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ajack/research/codes/examples/step-11

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ajack/research/codes/examples/step-11

# Include any dependencies generated for this target.
include CMakeFiles/step-11.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/step-11.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/step-11.dir/flags.make

CMakeFiles/step-11.dir/step-11.cc.o: CMakeFiles/step-11.dir/flags.make
CMakeFiles/step-11.dir/step-11.cc.o: step-11.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/ajack/research/codes/examples/step-11/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/step-11.dir/step-11.cc.o"
	/opt/apps/gcc/4.7.2/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/step-11.dir/step-11.cc.o -c /home/ajack/research/codes/examples/step-11/step-11.cc

CMakeFiles/step-11.dir/step-11.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/step-11.dir/step-11.cc.i"
	/opt/apps/gcc/4.7.2/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/ajack/research/codes/examples/step-11/step-11.cc > CMakeFiles/step-11.dir/step-11.cc.i

CMakeFiles/step-11.dir/step-11.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/step-11.dir/step-11.cc.s"
	/opt/apps/gcc/4.7.2/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/ajack/research/codes/examples/step-11/step-11.cc -o CMakeFiles/step-11.dir/step-11.cc.s

CMakeFiles/step-11.dir/step-11.cc.o.requires:
.PHONY : CMakeFiles/step-11.dir/step-11.cc.o.requires

CMakeFiles/step-11.dir/step-11.cc.o.provides: CMakeFiles/step-11.dir/step-11.cc.o.requires
	$(MAKE) -f CMakeFiles/step-11.dir/build.make CMakeFiles/step-11.dir/step-11.cc.o.provides.build
.PHONY : CMakeFiles/step-11.dir/step-11.cc.o.provides

CMakeFiles/step-11.dir/step-11.cc.o.provides.build: CMakeFiles/step-11.dir/step-11.cc.o

# Object files for target step-11
step__11_OBJECTS = \
"CMakeFiles/step-11.dir/step-11.cc.o"

# External object files for target step-11
step__11_EXTERNAL_OBJECTS =

step-11: CMakeFiles/step-11.dir/step-11.cc.o
step-11: CMakeFiles/step-11.dir/build.make
step-11: /opt/apps/gcc4_7/openmpi1_6/dealii/8.5.0/lib/libdeal_II.g.so.8.5.0
step-11: /opt/apps/gcc4_7/openmpi1_6/p4est/1.1/DEBUG/lib/libp4est.so
step-11: /opt/apps/gcc4_7/openmpi1_6/p4est/1.1/DEBUG/lib/libsc.so
step-11: /usr/lib64/libbz2.so
step-11: /opt/apps/gcc4_7/openmpi/1.6.5/lib/libmpi_f90.so
step-11: /opt/apps/gcc4_7/openmpi/1.6.5/lib/libmpi_f77.so
step-11: /usr/lib64/libz.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libpiro.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstokhos_sacado.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstokhos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libmoochothyra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libmoocho.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/librythmos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/liblocathyra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/liblocaepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/liblocalapack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libloca.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libnoxepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libnoxlapack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libnox.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libphalanx.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libintrepid.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteko.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikosbelos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikosaztecoo.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikosamesos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikosml.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libstratimikosifpack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libanasazitpetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libModeLaplace.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libanasaziepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libanasazi.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libbelostpetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libbelosepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libbelos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libml.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libifpack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libpamgen_extras.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libpamgen.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libamesos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libaztecoo.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libisorropia.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/liboptipack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libthyratpetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libthyraepetraext.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libthyraepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libthyracore.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libepetraext.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtpetrarti.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtpetraext.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtpetrainout.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtpetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtriutils.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libglobipack.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libshards.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libzoltan.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libepetra.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libsacado.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libkokkosdisttsqr.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libkokkosnodetsqr.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libkokkoslinalg.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libkokkosnodeapi.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libkokkos.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/librtop.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libtpi.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteuchosremainder.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteuchosnumerics.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteuchoscomm.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteuchosparameterlist.so
step-11: /opt/apps/gcc4_7/openmpi1_6/trilinos/11.12.1/lib/libteuchoscore.so
step-11: /opt/apps/gcc4_7/atlas/3.10.2/lib/libtatlas.so
step-11: /opt/apps/gcc4_7/openmpi/1.6.5/lib/libmpi_cxx.so
step-11: /opt/apps/gcc4_7/openmpi1_6/petsc/3.4.4/real-double/lib/libpetsc.so
step-11: /opt/apps/gcc4_7/atlas/3.10.1/lib/liblapack.a
step-11: /opt/apps/gcc4_7/atlas/3.10.1/lib/libcblas.a
step-11: /opt/apps/gcc4_7/atlas/3.10.1/lib/libf77blas.a
step-11: /opt/apps/gcc4_7/atlas/3.10.1/lib/libatlas.a
step-11: /usr/lib64/libX11.so
step-11: /opt/apps/gcc4_7/openmpi1_6/parmetis/4.0.3/lib/libparmetis.so
step-11: /opt/apps/gcc4_7/metis/5.0.2/lib/libmetis.a
step-11: /opt/apps/gcc4_7/openmpi1_6/phdf5/1.8.8/lib/libhdf5_hl.so
step-11: /opt/apps/gcc4_7/openmpi1_6/phdf5/1.8.8/lib/libhdf5.so
step-11: /opt/apps/gcc4_7/atlas/3.10.1/lib/libtatlas.so
step-11: /opt/apps/gcc4_7/openmpi/1.6.5/lib/libmpi.so
step-11: /usr/lib64/librdmacm.so
step-11: /usr/lib64/libibverbs.so
step-11: /opt/mellanox/mxm/lib/libmxm.so
step-11: /opt/torque/4.2.9/lib/libtorque.so
step-11: /usr/lib64/libnuma.so
step-11: /usr/lib64/libutil.so
step-11: CMakeFiles/step-11.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable step-11"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/step-11.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/step-11.dir/build: step-11
.PHONY : CMakeFiles/step-11.dir/build

CMakeFiles/step-11.dir/requires: CMakeFiles/step-11.dir/step-11.cc.o.requires
.PHONY : CMakeFiles/step-11.dir/requires

CMakeFiles/step-11.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/step-11.dir/cmake_clean.cmake
.PHONY : CMakeFiles/step-11.dir/clean

CMakeFiles/step-11.dir/depend:
	cd /home/ajack/research/codes/examples/step-11 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ajack/research/codes/examples/step-11 /home/ajack/research/codes/examples/step-11 /home/ajack/research/codes/examples/step-11 /home/ajack/research/codes/examples/step-11 /home/ajack/research/codes/examples/step-11/CMakeFiles/step-11.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/step-11.dir/depend
