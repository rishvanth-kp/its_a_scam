#
#   Copyright (C) 2019 Rishvanth Prabakar
#
#   Authors: Rish Prabakar
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

PROGS = rna_seq_qc

CXX = g++
CXXFLAGS = -Wall
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

all: $(PROGS)

$(PROGS): 

rna_seq_qc: rna_seq_qc.cpp GenomicRegion.o BedReader.o \
              GtfReader.o GenomicArray.o PreprocessGff.o \
							SamReader.o \
              StepVector.hpp GenomicStepVector.hpp
	$(CXX) $(CXXFLAGS) $^ -o $@


%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $<

%:%.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm $(PROGS) *.o

.PHONY: clean
