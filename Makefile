
CXX             =g++
FC              =g++
LD              =g++
MKDEP           =makedepend

WRKDIR          =$(shell pwd)
OBJDIR          =$(WRKDIR)/obj
INCDIR          =$(WRKDIR)/inc
SRCDIR          =$(WRKDIR)/src
INCROOT         =$(shell root-config --incdir)

CFLAGS  = -Dextname -I/usr/local/include

#CPPFLAGS        =-I /c/cern/rootcyg/include -I $(INCROOT)/ -I $(WRKDIR) -I $(INCDIR)/ -std=c++11
CPPFLAGS        =                           -I $(INCROOT) -I $(WRKDIR) -I $(INCDIR)/ -std=c++17
CXXFLAGS        =-fPIC -g -O
MKDEPFLAGS      =-Y ${INCL} -m -w 110

DICTKDetSim = KDetSimDictUX

KDetSimSL = KDetSim.sl

HDRSKDetSim =$(INCDIR)/KPad.h $(INCDIR)/KDetector.h $(INCDIR)/KField.h $(INCDIR)/K3D.h $(INCDIR)/KPixel.h $(INCDIR)/KGeometry.h $(INCDIR)/KMaterial.h $(INCDIR)/KStruct.h $(INCDIR)/KStrip.h  $(INCDIR)/KMesh.h  $(INCDIR)/KImplant3D.h $(INCDIR)/KImplant2D.h

DICTHDRSKDetSim  = $(HDRSKDetSim) $(INCDIR)/KDetSim_LinkDef.h

SRCSKDetSim      =$(SRCDIR)/KPad.cxx $(SRCDIR)/KDetector.cxx $(SRCDIR)/KField.cxx $(SRCDIR)/K3D.cxx $(SRCDIR)/KPixel.cxx $(SRCDIR)/KGeometry.cxx $(SRCDIR)/KMaterial.cxx $(SRCDIR)/KStruct.cxx $(SRCDIR)/KStrip.cxx $(SRCDIR)/KMesh.cxx $(SRCDIR)/KImplant3D.cxx $(SRCDIR)/KImplant2D.cxx

OBJSKDetSim      =$(OBJDIR)/KPad.o $(OBJDIR)/KDetector.o $(OBJDIR)/KField.o $(OBJDIR)/K3D.o $(OBJDIR)/KPixel.o $(OBJDIR)/KGeometry.o $(OBJDIR)/KMaterial.o $(OBJDIR)/KStruct.o $(OBJDIR)/KStrip.o $(OBJDIR)/KImplant2D.o $(OBJDIR)/KImplant3D.o $(OBJDIR)/KMesh.o $(OBJDIR)/nrutil.o $(OBJDIR)/$(DICTKDetSim).o


ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
LIBSLIN       =$(ROOTGLIBS) -lpthread -lm -ldl
LIBSMAC       =$(ROOTGLIBS) -lpthread -lm -ldl


KDetSim:  ${OBJSKDetSim} $(OBJDIR)/$(DICTKDetSim).o
	${LD} -shared ${CXXFLAGS} -o lib/${KDetSimSL} ${OBJSKDetSim} ${LIBSLIN}

# KDetSim
$(OBJDIR)/$(DICTKDetSim).o:  ${DICTHDRSKDetSim}
		rootcint -f $(SRCDIR)/${DICTKDetSim}.cxx -c ${DICTHDRSKDetSim}
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $(SRCDIR)/${DICTKDetSim}.cxx -o $(OBJDIR)/${DICTKDetSim}.o

clean:
		@rm ${OBJSKDetSim} lib/${KDetSimSL} $(SRCDIR)/${DICTKDetSim}.*

$(OBJDIR)/%.o:  $(SRCDIR)/%.cxx
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
