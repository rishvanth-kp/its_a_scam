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
          alignment_gc_distribution \
					alignment_depth \
					mismatch_profile coverage_test \
					genomic_feature_count mismatch_profile_vaf \
					mismatch_jaccard \
					gene_count add_barcode_to_fastq_name \
					barcode_count bc_in_region bc_tss_matrix \
					bc_feature_matrix bc_frag_dist genome_stats \
          tss_enrichment sc_rnaseq_feature_matrix \
          get_clipped_bases nametag_to_samtag \
          bc_mismatch_matrix get_reference_base

CXX = g++
CXXFLAGS = -Wall
LDLIBS = -lhts
OPTFLAGS = -O2
DEBUGFLAGS = -g
PROFFLAGS = -pg

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef PROF
CXXFLAGS += $(PROFFLAGS)
endif

all: $(PROGS)

$(PROGS):

rna_seq_qc: rna_seq_qc.cpp GenomicRegion.o BedReader.o \
              GtfReader.o PreprocessGff.o \
							SamReader.o SamEntry.o \
							AlignedGenomicFeature.o \
							AlignmentMismatches.o \
              StepVector.hpp GenomicStepVector.hpp

alignment_quality_distribution: SamReader.o SamEntry.o

alignment_length_distribution: SamReader.o SamEntry.o

alignment_gc_distribution: SamReader.o SamEntry.o

alignment_depth: SamReader.o SamEntry.o GenomicRegion.o StepVector.hpp \
									GenomicStepVector.hpp

mismatch_profile: SamReader.o SamEntry.o GenomicRegion.o FeatureVector.hpp \
									StepVector.hpp GenomicStepVector.hpp

mismatch_profile_vaf: SamReader.o SamEntry.o GenomicRegion.o \
											AlignmentMismatches.o FeatureVector.hpp \
											StepVector.hpp GenomicStepVector.hpp

mismatch_jaccard: SamReader.o SamEntry.o GenomicRegion.o \
										AlignmentMismatches.o FeatureVector.hpp \
										StepVector.hpp GenomicStepVector.hpp

coverage_test: SamReader.o SamEntry.o GenomicRegion.o FeatureVector.hpp \
									StepVector.hpp GenomicStepVector.hpp

genomic_feature_count: SamReader.o SamEntry.o AlignedGenomicFeature.o \
												StepVector.hpp GenomicStepVector.hpp \
												GtfReader.o GenomicRegion.o

gene_count: SamReader.o SamEntry.o GtfReader.o GeneCount.o StepVector.hpp \
							GenomicStepVector.hpp GenomicRegion.o FeatureVector.hpp

add_barcode_to_fastq_name: FastqReader.o

barcode_count: SamReader.o SamEntry.o 

bc_in_region: SamReader.o SamEntry.o BedReader.o StepVector.hpp \
                GenomicStepVector.hpp GenomicRegion.o

bc_tss_matrix: SamReader.o SamEntry.o BedReader.o StepVector.hpp \
                FeatureVector.hpp GenomicStepVector.hpp GenomicRegion.o

bc_feature_matrix: SamReader.o SamEntry.o StepVector.hpp BedReader.o \
	 									GenomicStepVector.hpp GtfReader.o GenomicRegion.o \
										AlignedOverlapGenomicFeature.o

bc_frag_dist: SamReader.o SamEntry.o

genome_stats: GenomeReader.o

tss_enrichment: SamReader.o SamEntry.o BedReader.o StepVector.hpp \
                  FeatureVector.hpp GenomicStepVector.hpp GenomicRegion.o

sc_rnaseq_feature_matrix: SamReader.o SamEntry.o AlignedGenomicFeature.o \
												    StepVector.hpp GenomicStepVector.hpp \
												    GtfReader.o GenomicRegion.o

get_clipped_bases: SamReader.o SamEntry.o

nametag_to_samtag: SamReader.o SamEntry.o

bc_mismatch_matrix: SamReader.o SamEntry.o GenomicRegion.o \
										  AlignmentMismatches.o FeatureVector.hpp \
										  StepVector.hpp GenomicStepVector.hpp

get_reference_base: GenomeReader.o BedReader.o 



%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

%:%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm $(PROGS) *.o

.PHONY: clean
