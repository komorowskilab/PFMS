# 
#  Copyright (C) 2009-2010
#  MOE Key Laboratory of Bioinformatics and 
#  Bioinformatics Division, TNLIST
#  Department of Automation, Tsinghua University
#  
#  Authors: Xi Wang and Xuegong Zhang
# 
#  This file is part of SeqSite.
# 
#  SeqSite is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  SeqSite is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

PROGS = SeqSite

BIN_DIR = ~/bin
HEADERS = $(shell ls *.h)
SOURCES = $(shell ls *.cpp)
DIR = $(shell basename `pwd`)

#OBJS = $(subst .cpp,.o,$(filter-out $(addsuffix .cpp,$(PROGS)),$(SOURCES)))

CXX = g++
CFLAGS = -Wall -fPIC
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

#ifdef OPT
CFLAGS += $(OPTFLAGS)
#endif


all:	$(PROGS)

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

#	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

SeqSite:     seqsite.o lm.o reads.o global.o bindSig.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

install: all
	install -m 755 $(PROGS) $(BIN_DIR)

uninstall:
	rm $(BIN_DIR)/SeqSite

tar: clean
	@tar czvf ../$(DIR).tar.gz ../$(DIR)

clean:
	@echo "clean compiled files"
	@-rm -f $(PROGS) *.o *.so *.a
