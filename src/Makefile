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

PROGS = rna_seq_qc alignment_quality_distribution \
          alignment_length_distribution \
          alignment_gc_distribution

CXX = g++
CXXFLAGS = -Wall
INCLUDEDIRS = -I/home/rish/code/tools/htslib/htslib_install/include/
LIBDIRS = -L/home/rish/code/tools/htslib/htslib_install/lib/
LIBS = -lhts
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
							SamReader.o SamEntry.o \
              StepVector.hpp GenomicStepVector.hpp

alignment_quality_distribution: SamReader.o SamEntry.o

alignment_length_distribution: SamReader.o SamEntry.o

alignment_gc_distribution: SamReader.o SamEntry.o

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< $(INCLUDEDIRS)

%:%.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLUDEDIRS) $(LIBDIRS) $(LIBS)

clean:
	rm $(PROGS) *.o

.PHONY: clean
