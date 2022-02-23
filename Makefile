# make all: compiles all
# clean: removes all .o and .mod files, all ~files etc.
# cleanall: removes also .eo files


#COMMENTS BEGIN
#NOTE: SINCE /usr/include/hdf5.mod WAS COMPILED WITH GFORTRAN, IFORT NEEDS A LOCAL LIBRARY!!!!
#FURTHER, NEED TO EXPORT LIBRARY PATH THE LIBRARY PATH:
#ADD_LIB_PATH=/home/moon/levin/lib64
#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NEW_LIB_PATH"
#export LD_LIBRARY_PATH
#COMMENTS END

.SUFFIXES: .f90


#################USER INPUT: CHOOSE COMPILER AND DEBUG OPTION############

#for now in shell script
COMPILER = gfortran
#COMPILER = ifort
#COMPILER = mpif90

#DEBUG = true

###########################HDF5 ENVIRONMENT VARIABLE#####################

ifndef DIR_HDF5
   $(info DIR_HDF5  Not found)
   $(error Undefined variable DIR_HDF5: export DIR_HDF5=<location of hdf5 lib>)
endif


#########################################################################

#operating system
UNAME := $(shell uname)

#####################GFORTRAN COMPILER OPTIONS ETC#######################

ifeq ($(COMPILER),gfortran)
   F90 = gfortran
   LD = gfortran

   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                      \
                   $(DIR_HDF5)/lib/libhdf5_fortran.so   \
                   $(DIR_HDF5)/lib/libhdf5.so           \
                   $(DIR_HDF5)/lib/libhdf5_hl.so        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                         \
                   $(DIR_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(DIR_HDF5)/lib/libhdf5.dylib           \
                   $(DIR_HDF5)/lib/libhdf5_hl.dylib        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp
#
   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -Og -Wall -Wextra -Warray-temporaries -Wconversion -pedantic-errors -fcheck=all -fbacktrace
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -ffree-line-length-none -c -O3 $(OMP_FLAG) $(DEBUG_FLAGS)
   MFLAGS = -J

   CDEFINED = true

endif

#####################IFORT COMPILER OPTIONS ETC#########################

ifeq ($(COMPILER),ifort)

   F90 = ifort
   LD = ifort

   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                      \
                   $(DIR_HDF5)/lib/libhdf5_fortran.so   \
                   $(DIR_HDF5)/lib/libhdf5.so           \
                   $(DIR_HDF5)/lib/libhdf5_hl.so        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                         \
                   $(DIR_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(DIR_HDF5)/lib/libhdf5.dylib           \
                   $(DIR_HDF5)/lib/libhdf5_hl.dylib        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp

   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -traceback -fpe0 -mcmodel=medium -shared-intel -check all -debug all
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -w -c -O3 -mp -assume byterecl -r8 $(OMP_FLAG) $(DEBUG_FLAG)
   MFLAGS = -module

   CDEFINED = true

endif

#####################MPIF909 COMPILER OPTIONS ETC########################

ifeq ($(COMPILER),mpif90)
   F90 = mpif90
   LD = mpif90


   ifeq ($(UNAME), Linux)
      $(info OS detected: Linux)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                      \
                   $(DIR_HDF5)/lib/libhdf5_fortran.so   \
                   $(DIR_HDF5)/lib/libhdf5.so           \
                   $(DIR_HDF5)/lib/libhdf5_hl.so        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.so
   else ifeq ($(UNAME), Darwin)
      $(info OS detected: MAC)
      DIR_MOD_HDF5=$(DIR_HDF5)/include
      DIR_LIB_HDF5=$(DIR_HDF5)/lib                         \
                   $(DIR_HDF5)/lib/libhdf5_fortran.dylib   \
                   $(DIR_HDF5)/lib/libhdf5.dylib           \
                   $(DIR_HDF5)/lib/libhdf5_hl.dylib        \
                   $(DIR_HDF5)/lib/libhdf5hl_fortran.dylib
   else
      $(error No OS detected, please have a look into the Makefile)
   endif

   OMP_FLAG=-fopenmp
#
   ifeq ($(DEBUG),true)
      DEBUG_FLAGS = -Og -Wall -Wextra -Warray-temporaries -Wconversion -pedantic-errors -fcheck=all -fbacktrace
   else
      DEBUG_FLAGS =
   endif

   CFLAGS = -ffree-line-length-none -c -O3 $(OMP_FLAG) $(DEBUG_FLAGS)
   MFLAGS = -J

   CDEFINED = true

endif

#########################################################################

ifneq ($(CDEFINED),true)
   ifneq ($(MAKECMDGOALS),clean)
      $(error compiler $(COMPILER) not defined)
   endif
endif

########################################################################

DIR_SRC = src
DIR_OBJ = objects/objects
DIR_MOD = modules/modules

DIR_SRC_DIFF1D = src_diff1d
DIR_OBJ_DIFF1D = objects/objects_diff1d
DIR_MOD_DIFF1D = modules/modules_diff1d

DIR_SRC_SC1D = src_sc1d
DIR_OBJ_SC1D = objects/objects_sc1d
DIR_MOD_SC1D = modules/modules_sc1d

DIR_SRC_SC2D = src_sc2d
DIR_OBJ_SC2D = objects/objects_sc2d
DIR_MOD_SC2D = modules/modules_sc2d

DIR_SRC_SC3D = src_sc3d
DIR_OBJ_SC3D = objects/objects_sc3d
DIR_MOD_SC3D = modules/modules_sc3d

DIR_SRC_MODEL = src_model
DIR_OBJ_MODEL = objects/objects_model
DIR_MOD_MODEL = modules/modules_model

DIR_SRC_OPAL = src_opal
DIR_OBJ_OPAL = objects/objects_opal
DIR_MOD_OPAL = modules/modules_opal

DIR_SRC_SURFB = src_surfb
DIR_OBJ_SURFB = objects/objects_surfb
DIR_MOD_SURFB = modules/modules_surfb

DIR_SRC_SC3DAMRVAC = src_sc3damrvac
DIR_OBJ_SC3DAMRVAC = objects/objects_sc3damrvac
DIR_MOD_SC3DAMRVAC = modules/modules_sc3damrvac

#module files
OBJSM_TYPE = $(DIR_OBJ)/mod_type.o

OBJSM_DIFF1D = $(DIR_OBJ)/mod_type.o                      \
               $(DIR_OBJ)/mod_integ1d.o         \
               $(DIR_OBJ)/mod_interp1d.o         \
               $(DIR_OBJ)/mod_interp2d.o     \
               $(DIR_OBJ)/mod_hello.o         \
               $(DIR_OBJ)/mod_math.o         \
               $(DIR_OBJ_OPAL)/mod_opal.o   \
               $(DIR_OBJ_DIFF1D)/mod_diff1d.o

OBJSM_SC1D = $(DIR_OBJ)/mod_type.o                      \
             $(DIR_OBJ)/mod_integ1d.o         \
             $(DIR_OBJ)/mod_interp1d.o         \
             $(DIR_OBJ)/mod_interp2d.o         \
             $(DIR_OBJ)/mod_hello.o         \
             $(DIR_OBJ)/mod_math.o         \
             $(DIR_OBJ)/mod_ng_extrapol.o      \
             $(DIR_OBJ)/mod_sparse.o      \
             $(DIR_OBJ_OPAL)/mod_opal.o    \
             $(DIR_OBJ_SC1D)/mod_sc1d.o

OBJSM_SC2D = $(DIR_OBJ)/mod_type.o                      \
             $(DIR_OBJ)/mod_integ1d.o         \
             $(DIR_OBJ)/mod_interp1d.o         \
             $(DIR_OBJ)/mod_interp2d.o         \
             $(DIR_OBJ)/mod_hello.o         \
             $(DIR_OBJ)/mod_math.o         \
             $(DIR_OBJ)/mod_ng_extrapol.o      \
             $(DIR_OBJ)/mod_sparse.o      \
             $(DIR_OBJ_OPAL)/mod_opal.o    \
             $(DIR_OBJ_SC2D)/mod_sc2d.o

OBJSM_SC3D = $(DIR_OBJ)/mod_type.o                      \
             $(DIR_OBJ)/mod_integ1d.o         \
             $(DIR_OBJ)/mod_interp1d.o         \
             $(DIR_OBJ)/mod_interp2d.o         \
             $(DIR_OBJ)/mod_interp3d.o         \
             $(DIR_OBJ)/mod_sort.o \
             $(DIR_OBJ)/mod_grid.o \
             $(DIR_OBJ)/mod_hello.o         \
             $(DIR_OBJ)/mod_math.o         \
             $(DIR_OBJ)/mod_ng_extrapol.o      \
             $(DIR_OBJ)/mod_sparse.o      \
             $(DIR_OBJ_SC3D)/mod_io.o   \
             $(DIR_OBJ_SC3D)/mod_grid3d.o   \
             $(DIR_OBJ_SC3D)/mod_model3d.o  \
             $(DIR_OBJ_SC3D)/mod_angles.o   \
             $(DIR_OBJ_SC3D)/mod_frequencies.o   \
             $(DIR_OBJ_SC3D)/mod_conttrans3d.o   \
             $(DIR_OBJ_SC3D)/mod_bcondition3d.o   \
             $(DIR_OBJ_SC3D)/mod_benchmark3d.o   \
             $(DIR_OBJ_OPAL)/mod_opal.o                  \
             $(DIR_OBJ_SC3D)/mod_scont_new3d.o      \
             $(DIR_OBJ_SC3D)/mod_params3d.o

OBJSM_MODEL = $(DIR_OBJ)/mod_type.o \
              $(DIR_OBJ)/mod_integ1d.o         \
              $(DIR_OBJ)/mod_interp1d.o               \
              $(DIR_OBJ)/mod_interp2d.o                  \
              $(DIR_OBJ)/mod_interp3d.o               \
              $(DIR_OBJ_OPAL)/mod_opal.o    \
              $(DIR_OBJ)/mod_sort.o  \
              $(DIR_OBJ)/mod_grid.o  \
              $(DIR_OBJ)/mod_amrvac_reader.o \
              $(DIR_OBJ)/mod_hello.o         \
              $(DIR_OBJ)/mod_math.o         \
              $(DIR_OBJ_MODEL)/mod_model.o

OBJSM_OPAL = $(DIR_OBJ)/mod_type.o \
             $(DIR_OBJ)/mod_interp2d.o \
             $(DIR_OBJ)/mod_hello.o         \
             $(DIR_OBJ)/mod_math.o         \
             $(DIR_OBJ_OPAL)/mod_opal.o

OBJSM_SURFB = $(DIR_OBJ)/mod_type.o \
              $(DIR_OBJ)/mod_integ1d.o         \
              $(DIR_OBJ)/mod_interp1d.o               \
              $(DIR_OBJ)/mod_interp2d.o               \
              $(DIR_OBJ)/mod_interp3d.o               \
              $(DIR_OBJ)/mod_hello.o         \
              $(DIR_OBJ)/mod_math.o         \
              $(DIR_OBJ_OPAL)/mod_opal.o \
              $(DIR_OBJ_SURFB)/mod_surfb.o

########################################################################

OBJS1_DIFF1D = $(DIR_OBJ)/model_laws.o       \
               $(DIR_OBJ)/info_region.o             \
               $(DIR_OBJ)/eispackEV.o
OBJS2_DIFF1D = $(OBJS1_DIFF1D) $(OBJSM_DIFF1D) $(DIR_OBJ_DIFF1D)/diff1d.o

OBJS1_SC1D = $(DIR_OBJ)/model_laws.o       \
             $(DIR_OBJ)/eispackEV.o        \
             $(DIR_OBJ)/formal_sc.o        \
             $(DIR_OBJ_SC1D)/mint_sc1d.o        \
             $(DIR_OBJ_SC1D)/formal_sc1d.o      \
             $(DIR_OBJ)/info_region.o             \
             $(DIR_OBJ_SC1D)/scont_new1d.o
OBJS2_SC1D = $(OBJS1_SC1D) $(OBJSM_SC1D) $(DIR_OBJ_SC1D)/sc1d.o

OBJS1_SC2D = $(DIR_OBJ)/model_laws.o       \
             $(DIR_OBJ)/eispackEV.o        \
             $(DIR_OBJ)/lapackINV.o        \
             $(DIR_OBJ)/formal_sc.o        \
             $(DIR_OBJ_SC2D)/mint_sc2d.o        \
             $(DIR_OBJ_SC2D)/formal_sc2d.o      \
             $(DIR_OBJ_SC2D)/scont_new2d.o      \
             $(DIR_OBJ)/info_region.o             \
             $(DIR_OBJ_SC2D)/output_sc2d.o
OBJS2_SC2D = $(OBJS1_SC2D) $(OBJSM_SC2D) $(DIR_OBJ_SC2D)/sc2d.o

OBJS1_SC3D = $(DIR_OBJ)/model_laws.o       \
             $(DIR_OBJ)/eispackEV.o        \
             $(DIR_OBJ)/lapackINV.o        \
             $(DIR_OBJ)/lebedev.o          \
             $(DIR_OBJ)/info_region.o      \
             $(DIR_OBJ)/formal_sc.o        \
             $(DIR_OBJ_SC3D)/debug_sc3d.o
OBJS2_SC3D = $(OBJS1_SC3D) $(OBJSM_SC3D) $(DIR_OBJ_SC3D)/sc3d.o


#for creating atmospheric model as input to sc3d.eo
OBJS1_MODEL = $(DIR_OBJ)/eispackEV.o              \
              $(DIR_OBJ)/info_region.o                   \
              $(DIR_OBJ)/model_laws.o       \
              $(DIR_OBJ_MODEL)/model3d.o                \
              $(DIR_OBJ_MODEL)/model_output.o

OBJS2_MODEL = $(OBJS1_MODEL) $(OBJSM_MODEL) $(DIR_OBJ_MODEL)/model.o


#for calculating surface brightness
OBJS1_SURFB = $(DIR_OBJ)/eispackEV.o              \
              $(DIR_OBJ)/info_region.o                   \
              $(DIR_OBJ)/model_laws.o       \
              $(DIR_OBJ)/get_iline.o              \
              $(DIR_OBJ_SURFB)/formal_lc.o                \
              $(DIR_OBJ_SURFB)/surfb_io.o

OBJS2_SURFB = $(OBJS1_SURFB) $(OBJSM_SURFB) $(DIR_OBJ_SURFB)/surfb.o

########################################################################

all: diff1d sc1d sc2d sc3d model surfb
#
########################################################################
#
diff1d: $(OBJS2_DIFF1D)
	 $(LD) $(OBJS2_DIFF1D) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o diff1d.eo

sc1d: $(OBJS2_SC1D)
	 $(LD) $(OBJS2_SC1D) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o sc1d.eo

sc2d: $(OBJS2_SC2D)
	 $(LD) $(OBJS2_SC2D) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o sc2d.eo

sc3d: $(OBJS2_SC3D)
	 $(LD) $(OBJS2_SC3D) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o sc3d.eo

model: $(OBJS2_MODEL)
	 $(LD) $(OBJS2_MODEL) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o model.eo

surfb: $(OBJS2_SURFB)
	 $(LD) $(OBJS2_SURFB) $(OMP_FLAG) -L $(DIR_LIB_HDF5) -o surfb.eo
#
########################################################################
#
$(DIR_OBJ)/mod_hello.o: $(DIR_SRC)/mod_hello.f90 \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_hello.o

$(DIR_OBJ)/mod_interp1d.o: $(DIR_SRC)/mod_interp1d.f90 \
                           $(DIR_OBJ)/mod_math.o \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp1d.o

$(DIR_OBJ)/mod_interp2d.o: $(DIR_SRC)/mod_interp2d.f90  \
                           $(DIR_OBJ)/mod_interp1d.o \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp2d.o

$(DIR_OBJ)/mod_interp3d.o: $(DIR_SRC)/mod_interp3d.f90  \
                           $(DIR_OBJ)/mod_interp1d.o \
                           $(DIR_OBJ)/info_region.o \
                           $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_interp3d.o

$(DIR_OBJ)/info_region.o: $(DIR_SRC)/info_region.f90  $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/info_region.o

$(DIR_OBJ)/mod_integ1d.o: $(DIR_SRC)/mod_integ1d.f90 \
                      $(DIR_OBJ)/mod_interp2d.o \
                      $(DIR_OBJ)/mod_math.o \
                      $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_integ1d.o

$(DIR_OBJ)/eispackEV.o: $(DIR_SRC)/eispackEV.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/eispackEV.o

$(DIR_OBJ)/get_iline.o: $(DIR_SRC)/get_iline.f90 $(DIR_OBJ)/model_laws.o $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/get_iline.o

$(DIR_OBJ)/lebedev.o: $(DIR_SRC)/lebedev.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/lebedev.o

$(DIR_OBJ)/lapackINV.o: $(DIR_SRC)/lapackINV.f
		$(F90) $(CFLAGS) $< -o $(DIR_OBJ)/lapackINV.o

$(DIR_OBJ)/mod_math.o: $(DIR_SRC)/mod_math.f90 \
                   $(DIR_OBJ)/eispackEV.o \
                   $(DIR_OBJ)/mod_sort.o \
                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_math.o

$(DIR_OBJ)/mod_sort.o: $(DIR_SRC)/mod_sort.f90 \
                       $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD)  $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_sort.o

$(DIR_OBJ)/mod_grid.o: $(DIR_SRC)/mod_grid.f90 \
                       $(DIR_OBJ)/mod_sort.o \
                       $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD)  $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_grid.o

$(DIR_OBJ)/mod_ng_extrapol.o: $(DIR_SRC)/mod_ng_extrapol.f90 \
                              $(DIR_OBJ)/mod_math.o \
                              $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_ng_extrapol.o

$(DIR_OBJ)/formal_sc.o: $(DIR_SRC)/formal_sc.f90 \
                        $(DIR_OBJ)/mod_integ1d.o \
                        $(DIR_OBJ)/mod_interp1d.o \
                        $(DIR_OBJ)/mod_interp2d.o \
                        $(DIR_OBJ)/mod_math.o \
                        $(DIR_OBJ)/mod_type.o 
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/formal_sc.o

$(DIR_OBJ)/model_laws.o: $(DIR_SRC)/model_laws.f90 \
                         $(DIR_OBJ)/mod_interp1d.o \
                         $(DIR_OBJ)/mod_interp2d.o \
                         $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $< -o $(DIR_OBJ)/model_laws.o

$(DIR_OBJ)/mod_sparse.o: $(DIR_SRC)/mod_sparse.f90 \
                         $(DIR_OBJ)/mod_ng_extrapol.o \
                         $(DIR_OBJ)/mod_math.o \
                         $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_sparse.o

$(DIR_OBJ)/mod_type.o: $(DIR_SRC)/mod_type.f90
		$(F90) $(CFLAGS) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD) $< -o $(DIR_OBJ)/mod_type.o

$(DIR_OBJ)/mod_amrvac_reader.o: $(DIR_SRC)/mod_amrvac_reader.f90 \
                                $(DIR_OBJ)/mod_sort.o \
                                $(DIR_OBJ)/mod_grid.o \
                                $(DIR_OBJ)/mod_interp1d.o \
                                $(DIR_OBJ)/mod_interp3d.o \
                                $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ)/mod_amrvac_reader.o


########################################################################

$(DIR_OBJ_OPAL)/mod_opal.o: $(DIR_SRC_OPAL)/mod_opal.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_OPAL) -I $(DIR_MOD) $< -o $(DIR_OBJ_OPAL)/mod_opal.o

########################################################################

$(DIR_OBJ_MODEL)/mod_model.o: $(DIR_SRC_MODEL)/mod_model.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_MODEL) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODEL)/mod_model.o

#$(DIR_OBJ_MODEL)/model1d.o: $(DIR_SRC_MODEL)/model1d.f90 $(DIR_OBJ)/model_laws.o $(OBJSM_MODEL)
#		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODEL)/model1d.o

#$(DIR_OBJ_MODEL)/model2d.o: $(DIR_SRC_MODEL)/model2d.f90 $(DIR_OBJ)/model_laws.o $(DIR_OBJ)/mod_interp1d.o $(OBJSM_MODEL)
#		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model2d.o

$(DIR_OBJ_MODEL)/model3d.o: $(DIR_SRC_MODEL)/model3d.f90 \
                            $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model3d.o

$(DIR_OBJ_MODEL)/model_output.o: $(DIR_SRC_MODEL)/model_output.f90 $(OBJSM_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_MODEL)/model_output.o

$(DIR_OBJ_MODEL)/model.o: $(DIR_SRC_MODEL)/model.f90 $(OBJSM_MODEL) $(OBJS1_MODEL)
		$(F90) $(CFLAGS) -I $(DIR_MOD_MODEL) -I $(DIR_MOD) $< -o $(DIR_OBJ_MODEL)/model.o

########################################################################

$(DIR_OBJ_SURFB)/mod_surfb.o: $(DIR_SRC_SURFB)/mod_surfb.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SURFB) -I $(DIR_MOD) $< -o $(DIR_OBJ_SURFB)/mod_surfb.o

$(DIR_OBJ_SURFB)/surfb_io.o: $(DIR_SRC_SURFB)/surfb_io.f90 $(OBJSM_SURFB) $(DIR_OBJ)/mod_math.o
		$(F90) $(CFLAGS) -I $(DIR_MOD_SURFB) -I $(DIR_MOD_OPAL) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SURFB)/surfb_io.o

$(DIR_OBJ_SURFB)/formal_lc.o: $(DIR_SRC_SURFB)/formal_lc.f90 $(OBJSM_SURFB)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SURFB) -I $(DIR_MOD) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SURFB)/formal_lc.o

$(DIR_OBJ_SURFB)/surfb.o: $(DIR_SRC_SURFB)/surfb.f90 $(OBJSM_SURFB) $(OBJS1_SURFB)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SURFB) -I $(DIR_MOD) $< -o $(DIR_OBJ_SURFB)/surfb.o

########################################################################

$(DIR_OBJ_DIFF1D)/diff1d.o: $(DIR_SRC_DIFF1D)/diff1d.f90 $(OBJSM_DIFF1D) $(OBJS1_DIFF1D)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_DIFF1D) $< -o $(DIR_OBJ_DIFF1D)/diff1d.o

$(DIR_OBJ_DIFF1D)/mod_diff1d.o: $(DIR_SRC_DIFF1D)/mod_diff1d.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS)  $(MFLAGS) $(DIR_MOD_DIFF1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_DIFF1D)/mod_diff1d.o


########################################################################

$(DIR_OBJ_SC1D)/mint_sc1d.o: $(DIR_SRC_SC1D)/mint_sc1d.f90 $(DIR_OBJ_SC1D)/formal_sc1d.o $(DIR_OBJ_SC1D)/scont_new1d.o $(OBJSM_SC1D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC1D)/mint_sc1d.o

$(DIR_OBJ_SC1D)/formal_sc1d.o: $(DIR_SRC_SC1D)/formal_sc1d.f90 $(DIR_OBJ)/mod_math.o $(OBJSM_SC1D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC1D)/formal_sc1d.o

$(DIR_OBJ_SC1D)/scont_new1d.o: $(DIR_SRC_SC1D)/scont_new1d.f90 \
                               $(DIR_OBJ)/mod_math.o \
                               $(DIR_OBJ)/mod_sparse.o $(OBJSM_SC1D) $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC1D)/scont_new1d.o

$(DIR_OBJ_SC1D)/sc1d.o: $(DIR_SRC_SC1D)/sc1d.f90 $(OBJSM_SC1D) $(OBJS1_SC1D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC1D)/sc1d.o

$(DIR_OBJ_SC1D)/mod_sc1d.o: $(DIR_SRC_SC1D)/mod_sc1d.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SC1D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC1D)/mod_sc1d.o

########################################################################

$(DIR_OBJ_SC2D)/mint_sc2d.o: $(DIR_SRC_SC2D)/mint_sc2d.f90 $(DIR_OBJ_SC2D)/formal_sc2d.o $(DIR_OBJ_SC2D)/scont_new2d.o $(OBJSM_SC2D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC2D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC2D)/mint_sc2d.o

$(DIR_OBJ_SC2D)/formal_sc2d.o: $(DIR_SRC_SC2D)/formal_sc2d.f90 $(DIR_OBJ)/mod_math.o $(OBJSM_SC2D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC2D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC2D)/formal_sc2d.o

$(DIR_OBJ_SC2D)/scont_new2d.o: $(DIR_SRC_SC2D)/scont_new2d.f90 $(DIR_OBJ)/mod_math.o $(DIR_OBJ)/mod_sparse.o $(OBJSM_SC2D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC2D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC2D)/scont_new2d.o

$(DIR_OBJ_SC2D)/sc2d.o: $(DIR_SRC_SC2D)/sc2d.f90 $(OBJSM_SC2D) $(OBJS1_SC2D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC2D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC2D)/sc2d.o

$(DIR_OBJ_SC2D)/mod_sc2d.o: $(DIR_SRC_SC2D)/mod_sc2d.f90 $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SC2D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC2D)/mod_sc2d.o

$(DIR_OBJ_SC2D)/output_sc2d.o: $(DIR_SRC_SC2D)/output_sc2d.f90 $(OBJSM_SC2D) 
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SC2D) -I $(DIR_MOD_HDF5) $< -o $(DIR_OBJ_SC2D)/output_sc2d.o


########################################################################

$(DIR_OBJ_SC3D)/mod_angles.o: $(DIR_SRC_SC3D)/mod_angles.f90 \
                              $(DIR_OBJ)/mod_math.o \
                              $(DIR_OBJ)/lebedev.o \
                              $(DIR_OBJ)/mod_integ1d.o \
                              $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_angles.o

$(DIR_OBJ_SC3D)/mod_frequencies.o: $(DIR_SRC_SC3D)/mod_frequencies.f90 \
                                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_frequencies.o

$(DIR_OBJ_SC3D)/mod_benchmark3d.o: $(DIR_SRC_SC3D)/mod_benchmark3d.f90 \
                                   $(DIR_OBJ_SC3D)/mod_conttrans3d.o \
                                   $(DIR_OBJ_SC3D)/mod_model3d.o \
                                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_benchmark3d.o

$(DIR_OBJ_SC3D)/debug_sc3d.o: $(DIR_SRC_SC3D)/debug_sc3d.f90 \
                              $(DIR_OBJ_SC3D)/mod_conttrans3d.o \
                              $(DIR_OBJ_SC3D)/mod_model3d.o \
                              $(OBJSM_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/debug_sc3d.o

$(DIR_OBJ_SC3D)/mod_scont_new3d.o: $(DIR_SRC_SC3D)/mod_scont_new3d.f90 \
                                   $(DIR_OBJ)/mod_math.o \
                                   $(DIR_OBJ)/mod_sparse.o \
                                   $(DIR_OBJ)/mod_ng_extrapol.o \
                                   $(DIR_OBJ_SC3D)/mod_params3d.o \
                                   $(DIR_OBJ_SC3D)/mod_grid3d.o \
                                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_scont_new3d.o

$(DIR_OBJ_SC3D)/sc3d.o: $(DIR_SRC_SC3D)/sc3d.f90 $(OBJSM_SC3D) $(OBJS1_SC3D)
		$(F90) $(CFLAGS) -I $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/sc3d.o

$(DIR_OBJ_SC3D)/mod_params3d.o: $(DIR_SRC_SC3D)/mod_params3d.f90 \
                            $(OBJSM_TYPE)
		$(F90) $(CFLAGS) $(MFLAGS) $(DIR_MOD_SC3D) -I $(DIR_MOD) $< -o $(DIR_OBJ_SC3D)/mod_params3d.o

$(DIR_OBJ_SC3D)/mod_bcondition3d.o: $(DIR_SRC_SC3D)/mod_bcondition3d.f90 \
                                    $(DIR_OBJ_SC3D)/mod_grid3d.o \
                                  $(OBJSM_TYPE)
		$(F90) $(CFLAGS)  -I $(DIR_MOD) -I $(DIR_MOD_SC3D) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_bcondition3d.o

$(DIR_OBJ_SC3D)/mod_model3d.o: $(DIR_SRC_SC3D)/mod_model3d.f90 \
                               $(DIR_OBJ)/model_laws.o \
                               $(DIR_OBJ)/mod_interp3d.o \
                               $(DIR_OBJ_SC3D)/mod_frequencies.o \
                               $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SC3D) -I $(DIR_MOD_OPAL) -I $(DIR_MOD_HDF5) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_model3d.o

$(DIR_OBJ_SC3D)/mod_grid3d.o: $(DIR_SRC_SC3D)/mod_grid3d.f90 \
                              $(DIR_OBJ)/mod_integ1d.o \
                              $(DIR_OBJ)/mod_grid.o \
                              $(DIR_OBJ_SC3D)/mod_params3d.o \
                              $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SC3D) -I $(DIR_MOD_HDF5) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_grid3d.o

$(DIR_OBJ_SC3D)/mod_io.o: $(DIR_SRC_SC3D)/mod_io.f90 \
                          $(DIR_OBJ_SC3D)/mod_params3d.o \
                          $(DIR_OBJ_SC3D)/mod_frequencies.o \
                          $(DIR_OBJ_SC3D)/mod_grid3d.o \
                          $(DIR_OBJ_SC3D)/mod_angles.o \
                          $(DIR_OBJ_SC3D)/mod_model3d.o \
                          $(DIR_OBJ_SC3D)/mod_conttrans3d.o \
                          $(DIR_OBJ_SC3D)/mod_bcondition3d.o \
                          $(DIR_OBJ_OPAL)/mod_opal.o \
                          $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SC3D) -I $(DIR_MOD_HDF5) -I $(DIR_MOD_OPAL) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_io.o

$(DIR_OBJ_SC3D)/mod_conttrans3d.o: $(DIR_SRC_SC3D)/mod_conttrans3d.f90 \
                                   $(DIR_OBJ_SC3D)/mod_params3d.o \
                                   $(DIR_OBJ_SC3D)/mod_frequencies.o \
                                   $(DIR_OBJ_SC3D)/mod_grid3d.o \
                                   $(DIR_OBJ_SC3D)/mod_bcondition3d.o \
                                   $(DIR_OBJ_SC3D)/mod_angles.o \
                                   $(DIR_OBJ_SC3D)/mod_model3d.o \
                                   $(DIR_OBJ_SC3D)/mod_scont_new3d.o \
                                   $(DIR_OBJ)/mod_interp2d.o \
                                   $(DIR_OBJ)/mod_math.o \
                                   $(OBJSM_TYPE)
		$(F90) $(CFLAGS) -I $(DIR_MOD) -I $(DIR_MOD_SC3D) -I $(DIR_MOD_HDF5) -I $(DIR_MOD_OPAL) $(MFLAGS) $(DIR_MOD_SC3D) $< -o $(DIR_OBJ_SC3D)/mod_conttrans3d.o


########################################################################

# inputparam.mod gridstratrev.mod comvelo.mod comtluc.mod
.PHONY: clean cleanall

clean:	
		rm -f *.o *~ *.mod \#*
		rm -f $(DIR_OBJ)/*.o $(DIR_OBJ)/*~ $(DIR_OBJ)/*.mod
		rm -f $(DIR_MOD)/*.o $(DIR_MOD)/*~ $(DIR_MOD)/*.mod
		rm -f $(DIR_OBJ_DIFF1D)/*.o $(DIR_OBJ_DIFF1D)/*~ $(DIR_OBJ_DIFF1D)/*.mod
		rm -f $(DIR_MOD_DIFF1D)/*.o $(DIR_MOD_DIFF1D)/*~ $(DIR_MOD_DIFF1D)/*.mod
		rm -f $(DIR_OBJ_SC1D)/*.o $(DIR_OBJ_SC1D)/*~ $(DIR_OBJ_SC1D)/*.mod
		rm -f $(DIR_MOD_SC1D)/*.o $(DIR_MOD_SC1D)/*~ $(DIR_MOD_SC1D)/*.mod
		rm -f $(DIR_OBJ_SC2D)/*.o $(DIR_OBJ_SC2D)/*~ $(DIR_OBJ_SC2D)/*.mod
		rm -f $(DIR_MOD_SC2D)/*.o $(DIR_MOD_SC2D)/*~ $(DIR_MOD_SC2D)/*.mod
		rm -f $(DIR_OBJ_SC3D)/*.o $(DIR_OBJ_SC3D)/*~ $(DIR_OBJ_SC3D)/*.mod
		rm -f $(DIR_MOD_SC3D)/*.o $(DIR_MOD_SC3D)/*~ $(DIR_MOD_SC3D)/*.mod
		rm -f $(DIR_OBJ_MODEL)/*.o $(DIR_OBJ_MODEL)/*~ $(DIR_OBJ_MODEL)/*.mod
		rm -f $(DIR_MOD_MODEL)/*.o $(DIR_MOD_MODEL)/*~ $(DIR_MOD_MODEL)/*.mod
		rm -f $(DIR_OBJ_SURFB)/*.o $(DIR_OBJ_SURFB)/*~ $(DIR_OBJ_SURFB)/*.mod
		rm -f $(DIR_MOD_SURFB)/*.o $(DIR_MOD_SURFB)/*~ $(DIR_MOD_SURFB)/*.mod

		rm -f $(DIR_SRC)/*~
		rm -f $(DIR_SRC_DIFF1D)/*~
		rm -f $(DIR_SRC_SC1D)/*~
		rm -f $(DIR_SRC_SC2D)/*~
		rm -f $(DIR_SRC_SC3D)/*~
		rm -f $(DIR_SRC_MODEL)/*~

cleanall: 
		rm -f *.eo *.o *~ *.mod \#*

