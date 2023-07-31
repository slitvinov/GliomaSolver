.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

LIBS = \
	-fopenmp\
	-fpermissive\
	tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_release/libtbb.so\
	-ldl \
	-LmyVTK/lib/vtk-5.2/\
	-lm \
	-lpthread \
	-ltbb \
	-ltbbmalloc \
	-lvtkCommon \
	-lvtkDICOMParser \
	-lvtkexpat \
	-lvtkFiltering \
	-lvtkftgl \
	-lvtkGraphics \
	-lvtkHybrid \
	-lvtkImaging \
	-lvtkInfovis \
	-lvtkIO \
	-lvtkjpeg \
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

VPATH = Glioma/ MRAG/ MRAG/MRAGcore Glioma/Tests/

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
	-IGlioma/\
	-IGlioma/Operators/\
	-IGlioma/Tests/\
	-IMRAG/\
	-O3 \
	-Wno-deprecated\
	-ImyVTK/include/vtk-5.2/\
	-Itbb40_20120613oss/include/\

brain: $(OBJECTS); $(CXX) $^ -o $@ $(LIBS)
.cpp.o: $(CXX) $(CPPFLAGS) -c $^ -o $@
clean:; rm -rf $(OBJECTS) brain
