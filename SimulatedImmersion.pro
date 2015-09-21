CONFIG   -= app_bundle
#CONFIG -= qt

TARGET = SimulatedImmersion3D
# Directories
DESTDIR = dist/
OBJECTS_DIR = build/


QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3


SOURCES += \
    source/main.cpp \
    source/Mesh.cpp \
    source/Matrix.cpp \
    source/IO.cpp

HEADERS += \
    source/Vertex.h \
    source/Triangle.h \
    source/Tetra.h \
    source/Mesh.h \
    source/LibraryHeader.h \
    source/IO.h \
    source/Edge.h \
    source/Define.h \
    source/SuperTetras.h

