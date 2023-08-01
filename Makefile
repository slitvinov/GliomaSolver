.SUFFIXES:
.SUFFIXES: .cpp
.SUFFIXES: .o

LIBS_TBB = `pkg-config --libs tbb`
CXXFLAGS_TBB = `pkg-config --cflags tbb`

CXXFLAGS_VTK = -IVTK/include/vtk-5.2
LIBS_VTK = -LVTK/lib/vtk-5.2

LIBS = \
	$(LIBS_TBB)\
	$(LIBS_VTK)\
	-fopenmp\
	-fpermissive\
	-ltbbmalloc\
	-lvtkCommon \
	-lvtkDICOMParser\
	-lvtkexpat \
	-lvtkFiltering\
	-lvtkftgl\
	-lvtkIO\
	-lvtkjpeg\
	-lvtkmetaio\
	-lvtkNetCDF\
	-lvtkpng\
	-lvtksqlite\
	-lvtksys\
	-lvtktiff\
	-lvtkverdict\
	-lvtkzlib\

VPATH = Glioma MRAG MRAG/MRAGcore

OBJECTS = \
Glioma_main.o \
Glioma_ReactionDiffusion.o\
MRAGBoundaryBlockInfo.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \

CXXFLAGS = \
	-D_BLOCKSIZE_=16\
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
	-O3\
	-Wno-deprecated\
	$(CXXFLAGS_VTK)\
	$(CXXFLAGS_TBB)\

brain: $(OBJECTS); $(CXX) $(OBJECTS) -o $@ $(LIBS)
.cpp.o:; $(CXX) $(CXXFLAGS) -c $< -o $@
clean:; rm -rf $(OBJECTS) brain
