# Generated by CMake 2.6-patch 4

IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
  # Information for CMake 2.6 and above.
  SET("vtkCommon_LIB_DEPENDS" "general;vtksys;general;-lpthread;general;dl;general;-lm;")
  SET("vtkFiltering_LIB_DEPENDS" "general;vtkCommon;")
  SET("vtkGenericFiltering_LIB_DEPENDS" "general;vtkFiltering;general;vtkGraphics;")
  SET("vtkGraphics_LIB_DEPENDS" "general;vtkFiltering;general;vtkverdict;")
  SET("vtkHybrid_LIB_DEPENDS" "general;vtkRendering;general;vtkIO;general;vtkexoIIc;")
  SET("vtkIO_LIB_DEPENDS" "general;vtkFiltering;general;vtkDICOMParser;general;vtkNetCDF;general;vtkmetaio;general;vtksqlite;general;vtkpng;general;vtkzlib;general;vtkjpeg;general;vtktiff;general;vtkexpat;general;vtksys;")
  SET("vtkImaging_LIB_DEPENDS" "general;vtkFiltering;")
  SET("vtkInfovis_LIB_DEPENDS" "general;vtkWidgets;general;vtklibxml2;general;vtkalglib;")
  SET("vtkRendering_LIB_DEPENDS" "general;vtkGraphics;general;vtkImaging;general;vtkIO;general;vtkftgl;general;vtkfreetype;general;/usr/lib64/libGL.so;general;/usr/lib64/libXt.so;general;/usr/lib64/libSM.so;general;/usr/lib64/libICE.so;general;/usr/lib64/libX11.so;general;/usr/lib64/libXext.so;")
  SET("vtkViews_LIB_DEPENDS" "general;vtkInfovis;")
  SET("vtkVolumeRendering_LIB_DEPENDS" "general;vtkRendering;general;vtkIO;")
  SET("vtkWidgets_LIB_DEPENDS" "general;vtkRendering;general;vtkHybrid;")
  SET("vtkexoIIc_LIB_DEPENDS" "general;vtkNetCDF;")
  SET("vtkftgl_LIB_DEPENDS" "general;/usr/lib64/libGL.so;general;vtkfreetype;")
  SET("vtklibxml2_LIB_DEPENDS" "general;vtkzlib;general;dl;general;-lpthread;general;dl;general;m;")
  SET("vtkmetaio_LIB_DEPENDS" "general;vtkzlib;general;vtksys;")
  SET("vtkpng_LIB_DEPENDS" "general;vtkzlib;")
  SET("vtkproj4_LIB_DEPENDS" "general;m;")
  SET("vtksys_LIB_DEPENDS" "general;dl;")
  SET("vtktiff_LIB_DEPENDS" "general;vtkzlib;general;vtkjpeg;")
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
  # Information for CMake 2.4 and lower.
  SET("vtkCommon_LIB_DEPENDS" "vtksys;-lpthread;dl;-lm;")
  SET("vtkFiltering_LIB_DEPENDS" "vtkCommon;")
  SET("vtkGenericFiltering_LIB_DEPENDS" "vtkFiltering;vtkGraphics;")
  SET("vtkGraphics_LIB_DEPENDS" "vtkFiltering;vtkverdict;")
  SET("vtkHybrid_LIB_DEPENDS" "vtkRendering;vtkIO;vtkexoIIc;")
  SET("vtkIO_LIB_DEPENDS" "vtkFiltering;vtkDICOMParser;vtkNetCDF;vtkmetaio;vtksqlite;vtkpng;vtkzlib;vtkjpeg;vtktiff;vtkexpat;vtksys;")
  SET("vtkImaging_LIB_DEPENDS" "vtkFiltering;")
  SET("vtkInfovis_LIB_DEPENDS" "vtkWidgets;vtklibxml2;vtkalglib;")
  SET("vtkRendering_LIB_DEPENDS" "vtkGraphics;vtkImaging;vtkIO;vtkftgl;vtkfreetype;/usr/lib64/libGL.so;/usr/lib64/libXt.so;/usr/lib64/libSM.so;/usr/lib64/libICE.so;/usr/lib64/libX11.so;/usr/lib64/libXext.so;")
  SET("vtkViews_LIB_DEPENDS" "vtkInfovis;")
  SET("vtkVolumeRendering_LIB_DEPENDS" "vtkRendering;vtkIO;")
  SET("vtkWidgets_LIB_DEPENDS" "vtkRendering;vtkHybrid;")
  SET("vtkexoIIc_LIB_DEPENDS" "vtkNetCDF;")
  SET("vtkftgl_LIB_DEPENDS" "/usr/lib64/libGL.so;vtkfreetype;")
  SET("vtklibxml2_LIB_DEPENDS" "vtkzlib;dl;-lpthread;dl;m;")
  SET("vtkmetaio_LIB_DEPENDS" "vtkzlib;vtksys;")
  SET("vtkpng_LIB_DEPENDS" "vtkzlib;")
  SET("vtkproj4_LIB_DEPENDS" "m;")
  SET("vtksys_LIB_DEPENDS" "dl;")
  SET("vtktiff_LIB_DEPENDS" "vtkzlib;vtkjpeg;")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.4)
