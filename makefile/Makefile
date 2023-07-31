.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

LIBS = \
	-fopenmp\
	-fpermissive\
	../INSTALL/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/libtbb.so\
	-ldl \
	-L../INSTALL/myVTK/lib/vtk-5.2/\
	-lm \
	-lpthread \
	-ltbb \
	-ltbbmalloc \
	-lvtkCommon \
	-lvtkDICOMParser \
	-lvtkexoIIc \
	-lvtkexpat \
	-lvtkFiltering \
	-lvtkfreetype \
	-lvtkftgl \
	-lvtkGraphics \
	-lvtkHybrid \
	-lvtkImaging \
	-lvtkInfovis \
	-lvtkIO \
	-lvtkjpeg \
	-lvtklibxml2 \
	-lvtkmetaio \
	-lvtkNetCDF \
	-lvtkpng \
	-lvtkRendering \
	-lvtksqlite \
	-lvtksys \
	-lvtktiff \
	-lvtkverdict \
	-lvtkViews \
	-lvtkWidgets \
	-lvtkzlib \

VPATH = ../Glioma/ ../MRAG/ ../MRAG/MRAGcore ../Glioma/Tests/

OBJECTS = \
dat2VP.o \
Glioma_ComputePFF_CahnHilliard.o\
Glioma_main.o \
Glioma_ReactionDiffusion.o\
Glioma_UQ_DataPreprocessing.o\
MRAGBoundaryBlockInfo.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \
Test.o\

CPPFLAGS = \
	-D_BLOCKSIZE_=8\
	-D_BPD_=16\
	-D_DIM=3\
	-D_FMMSILENT\
	-D _MAXLEVEL_=4\
	-DNDEBUG\
	-D_RESJUMP_=1\
	-I../Glioma/\
	-I../Glioma/Operators/\
	-I../Glioma/Tests/\
	-I../MRAG/\
	-O3 \
	-Wno-deprecated\
	-I../INSTALL/myVTK/include/vtk-5.4/\
	-I../INSTALL/tbb40_20120613oss/include/\

brain: $(OBJECTS); $(CXX) $^ -o $@ $(LIBS)
.cpp.o: $(CXX) $(CPPFLAGS) -c $^ -o $@
