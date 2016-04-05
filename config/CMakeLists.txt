#download test data
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/data/NA12878.chrom20.bam")
	   MESSAGE("Downloading BAM")

	   file(DOWNLOAD
			ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
			${PROJECT_SOURCE_DIR}/data/NA12878.chrom20.bam
			SHOW_PROGRESS)
endif()
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/data/NA12878.chrom20.bam.bai")
	   MESSAGE("Downloading BAM Index")
	   file(DOWNLOAD
			ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.bai
			${PROJECT_SOURCE_DIR}/data/NA12878.chrom20.bam.bai
			SHOW_PROGRESS)
endif()

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/DataConfig.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/DataConfig.h"
)
