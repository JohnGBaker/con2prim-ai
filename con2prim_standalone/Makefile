# Set C++ compiler
CXX      = /opt/local/bin/g++-mp-9

# Set C++ compiler flags
CXXFLAGS = -Wall -fopenmp -O2 -march=native -DENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA

# Source file(s)
SRC = driver_conserv_to_prims.C

# Include file(s) (Note: in IGM most files are included in
# driver files and therefore there will be .c and .C files
# in the include list)
INC = standalone_conserv_to_prims_main_function_random_data.h \
      IllinoisGRMHD_headers.h \
      harm_primitives_headers.h \
      harm_u2p_util.c \
      harm_utoprim_2d.c \
      harm_utoprim_1d.c \
      harm_utoprim_1d_ee.c \
      harm_utoprim_1d_ee2.c \
      IllinoisGRMHD_EoS_lowlevel_functs.C \
      apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C \
      convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij.C \
      eigen.C \
      font_fix_hybrid_EOS.C \
      harm_primitives_lowlevel.C \
      inlined_functions.C

# Set executable name
EXECUTABLE = con2prim_random_input

# This will keep the Makefile working even if files "all" and "clean" exist
.PHONY: all clean

# All: default target
all: $(EXECUTABLE)

# Compile the standalone
$(EXECUTABLE): $(SRC) $(INC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

# Cleanup
clean:
	rm -rf $(EXECUTABLE)
