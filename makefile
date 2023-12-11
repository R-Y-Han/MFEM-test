# Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# Use the MFEM build directory
MFEM_DIR ?= /home/han/Code/MFEM/mfem
MFEM_BUILD_DIR ?= /home/han/Code/MFEM/mfem
SRC = ./
# dir for *.o
OBJ = ./obj/
# dir for lib
LIB = ./lib/
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_INSTALL_DIR = ../mfem
# CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

SEQ_EXAMPLES = HeatConduction

ifeq ($(MFEM_USE_MPI),NO)
   EXAMPLES = $(SEQ_EXAMPLES)
else
   EXAMPLES = $(PAR_EXAMPLES) $(SEQ_EXAMPLES)
endif
SUBDIRS =


SUBDIRS_ALL = $(addsuffix /all,$(SUBDIRS))
SUBDIRS_TEST = $(addsuffix /test,$(SUBDIRS))
SUBDIRS_CLEAN = $(addsuffix /clean,$(SUBDIRS))
SUBDIRS_TPRINT = $(addsuffix /test-print,$(SUBDIRS))

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec lib

# Remove built-in rule
%: %.cpp

# Replace the default implicit rule for *.cpp files
%: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

# Replace the default implicit rule for *.hpp files
$(OBJ)%.o: $(SRC)%.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(MFEM_FLAGS) $< -o $@ $(MFEM_LIBS)

# define src and obj for lib files
LIB_SRCS = $(wildcard $(SRC)*.hpp)
LIB_OBJS = $(patsubst $(SRC)%.hpp,$(OBJ)%.o,$(LIB_SRCS))
$(LIB)MY_MFEM_LIB.a: $(LIB_OBJS)
	ar rcs $@ $^

all: $(EXAMPLES) $(SUBDIRS_ALL)
lib: $(LIB)MY_MFEM_LIB.a

.PHONY: $(SUBDIRS_ALL) $(SUBDIRS_TEST) $(SUBDIRS_CLEAN) $(SUBDIRS_TPRINT)
$(SUBDIRS_ALL) $(SUBDIRS_TEST) $(SUBDIRS_CLEAN):
	$(MAKE) -C $(@D) $(@F)
$(SUBDIRS_TPRINT):
	@$(MAKE) -C $(@D) $(@F)

# Additional dependencies
HeatConduction: $(SRC)HeatConduction.hpp


MFEM_TESTS = EXAMPLES
include $(MFEM_TEST_MK)
test: $(SUBDIRS_TEST)
test-print: $(SUBDIRS_TPRINT)

# Testing: Parallel vs. serial runs
RUN_MPI = $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(MFEM_MPI_NP)
%-test-par: %
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example)
%-test-seq: %
	@$(call mfem-test,$<,, Serial example)
%-test-par-cuda: %
	@$(call mfem-test,$<, $(RUN_MPI), Parallel CUDA example,-d cuda)
%-test-seq-cuda: %
	@$(call mfem-test,$<,, Serial CUDA example,-d cuda)
%-test-par-hip: %
	@$(call mfem-test,$<, $(RUN_MPI), Parallel HIP example,-d hip)
%-test-seq-hip: %
	@$(call mfem-test,$<,, Serial HIP example,-d hip)

# Testing: Specific execution options
ex0-test-seq: ex0
	@$(call mfem-test,$<,, Serial example,,1)
ex0p-test-par: ex0p
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example,,1)
ex1-test-seq: ex1
	@$(call mfem-test,$<,, Serial example)
ex1p-test-par: ex1p
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example)
ex10-test-seq: ex10
	@$(call mfem-test,$<,, Serial example,-tf 5)
ex10p-test-par: ex10p
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example,-tf 5)
ex15-test-seq: ex15
	@$(call mfem-test,$<,, Serial example,-e 1)
ex15p-test-par: ex15p
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example,-e 1)
ex27-test-seq: ex27
	@$(call mfem-test,$<,, Serial example,-dg)
ex27p-test-par: ex27p
	@$(call mfem-test,$<, $(RUN_MPI), Parallel example,-dg)
# Testing: optional tests
ifeq ($(MFEM_USE_STRUMPACK),YES)
ex11p-test-strumpack: ex11p
	@$(call mfem-test,$<, $(RUN_MPI), STRUMPACK example,--strumpack)
test-par-YES: ex11p-test-strumpack
endif
ifeq ($(MFEM_USE_SUPERLU),YES)
ex11p-test-superlu: ex11p
	@$(call mfem-test,$<, $(RUN_MPI), SuperLU_DIST example,--superlu)
test-par-YES: ex11p-test-superlu
endif
ifeq ($(MFEM_USE_MKL_CPARDISO),YES)
ex11p-test-cpardiso: ex11p
	@$(call mfem-test,$<, $(RUN_MPI), MKL_CPARDISO example,--cpardiso)
test-par-YES: ex11p-test-cpardiso
endif

# Testing: "test" target and mfem-test* variables are defined in config/test.mk

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec $(SUBDIRS_CLEAN)

clean-build:
	rm -f *.o *~ $(SEQ_EXAMPLES) $(PAR_EXAMPLES)
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -f refined.mesh displaced.mesh mesh.* ex5.mesh ex6p-checkpoint.*
	@rm -rf Example5* Example9* Example15* Example16* Example23* ParaView
	@rm -f sphere_refined.* sol.* sol_u.* sol_p.* sol_r.* sol_i.*
	@rm -f ex9.mesh ex9-mesh.* ex9-init.* ex9-final.*
	@rm -f deformed.* velocity.* elastic_energy.* mode_* mode_deriv_* flux.*
	@rm -f ex5-p-*.bp ex9-p-*.bp ex12-p-*.bp ex16-p-*.bp
	@rm -f ex16.mesh ex16-mesh.* ex16-init.* ex16-final.*
	@rm -f vortex-mesh.* vortex.mesh vortex-?-init.* vortex-?-final.*
	@rm -f deformation.* pressure.*
	@rm -f ex20.dat ex20p_?????.dat gnuplot_ex20.inp gnuplot_ex20p.inp
	@rm -f ex21*.mesh ex21*.sol ex21p_*.*
	@rm -f ex23.mesh ex23-*.gf
	@rm -f ex25.mesh ex25-*.gf ex25p-*.*
	@rm -rf ex28_* ex28p_*
	@rm -rf cond.* cond_mesh.* cond_j.* dsol.* port_mesh.* port_mode.*