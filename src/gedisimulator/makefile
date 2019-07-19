# Makefile for GEDI simulator tools
LIBS = -lm -lgdal -lgeotiff 
BIN=
ARCH=$(shell uname -m)
LINKER=
COPY=cp
THIS=gediRat
OUTPUT=/usr/bin/$(THIS)

ifeq ("$(OS)","Windows_NT")
	LIBS = -lgdal_i -lgeotiff_i
	CFLAGS += -D_USE_MATH_DEFINES=1 -D_WIN32 -DWIN32 -DH5_BUILT_AS_DYNAMIC_LIB -DGSL_DLL
	BIN=.exe
	ARCH=$(PROCESSOR_ARCHITECTURE)
	LINKER=$(where link)
	COPY=copy /Y
	OUTPUT = $(CONDA_PREFIX)\Library\bin\$(THIS)$(BIN)
endif


LIBS += -lgsl -lgslcblas -ltiff -lhdf5 -L${GSL_ROOT} -L${HDF5_LIB}/lib 
CFLAGS += -I/usr/local/include -I$(HANCOCKTOOLS_ROOT) -I$(CMPFIT_ROOT) -I${LIBCLIDAR_ROOT} -I${HDF5_LIB}/include -I. -I/usr/include/libgeotiff -I/usr/include/gdal
CFLAGS += -Wall 
#CFLAGS += -Wl,--verbose
CFLAGS += -O3
#CFLAGS += -g
src = $(LIBCLIDAR_ROOT)/libLasProcess.c $(LIBCLIDAR_ROOT)/libLasRead.c $(LIBCLIDAR_ROOT)/tiffWrite.c $(LIBCLIDAR_ROOT)/gaussFit.c $(LIBCLIDAR_ROOT)/libLidVoxel.c  $(LIBCLIDAR_ROOT)/libTLSread.c  $(LIBCLIDAR_ROOT)/libLidarHDF.c $(GEDIRAT_ROOT)/gediIO.c $(LIBCLIDAR_ROOT)/libOctree.c $(GEDIRAT_ROOT)/gediNoise.c #$(GEDIRAT_ROOT)/photonCount.c
src += $(CMPFIT_ROOT)/mpfit.c
obj = $(src:.c=.o)
EXE=$(THIS)$(BIN)

# GSLFit=linear.o
# MIN=mpfit
# MPFIT=${CMPFIT_ROOT}/$(MIN)

CC=gcc

$(EXE): $(THIS).o $(obj)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) $(INCLS)

clean:
	rm -f *% *~ *.o

install: $(THIS)$(BIN)
	$(COPY) "$(THIS)$(BIN)" "$(OUTPUT)"

.PHONY: clean install