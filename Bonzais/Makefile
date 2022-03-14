# Squeezer Makefile
# Author: Ph. Gras CEA/IRFU Saclay
# Jan. 3, 11

.PHONY: all clean

SOURCES=pruner.cc Pruner.cc ShearsTChain.cc
OBJS=$(addsuffix .o, $(basename $(SOURCES)))
DEPS=$(addsuffix .d, $(basename $(SOURCES)))
SO=$(addsuffix .so, $(shell cat pruners.txt))

CPPFLAGS+=-MD -Wall
CXXFLAGS+=$(shell root-config --cflags) -I.
CXXFLAGS+=-g -fPIC
LDFLAGS+=$(shell root-config --libs --ldflags)
LINK.o = $(CXX) $(LDFLAGS) $(TARGET_ARCH)
ROOTCINT=rootcint
CINT_INCDIR=$(shell root-config --incdir)
EXE=pruner
all: $(EXE) $(SO)

clean:
	$(RM) $(EXE) *.o *.d *.so

$(EXE): $(OBJS)
	$(LINK.cc) $^ $(LOADLIBES) $(LDLIBS) -o $@

%Dict.o: %Dict.cc

%Dict.cc : %.h  %LinkDef.h
	$(ROOTCINT) -f $@ -c -I. -I$(CINT_INCDIR) $(DEFINES) $^

%.so: %.o
	$(LINK.cc) $(LOADLIBES) $(LDLIBS) -shared -o $@ $^

echo_%:
	@echo $* = $($*)

-include $(DEPS)

#ignore missing header file to prevent 'make' error after a header file
 #removal:
%.h:;
