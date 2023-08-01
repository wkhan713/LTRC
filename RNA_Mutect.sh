import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:CramToBam.wdl_copy/versions/4/plain-WDL/descriptor' as CramToBam
import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:BamRealigner_wdl_v1_0/versions/8/plain-WDL/descriptor' as BamRealigner
import 'https://api.firecloud.org/ga4gh/v1/tools/bbecerra:RemoveAtypicalReads/versions/8/plain-WDL/descriptor' as RemoveAtypicalReads

#The next 3 tasks correspond to the CallSomaticMutations_V131 task under the "Call Somatic Mutations For Capture" workflow.
task CallSomaticMutations_131_Prepare_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible
        Array[Int] continueOnReturnCode

		#TASK INPUT PARAMS
		File? sampleTargetsIntervalList
		File workspaceTargetsIntervalList
		File refFasta
		File refFastaDict
		File refFastaIdx
		Int nWay
		Float tBamSize
		Float tBaiSize
		Float nBamSize
		Float nBaiSize
		Float dbSNPVCFSize
		Float cosmicVCFSize
		Float normalPanelSize
		Float fastaSize
		Float fastaDictSize
	command <<<
		#increase verbosity
		set -x

		#Calculate disk size for all shards of the mutect run
		SIZE_FILE=split_base_sizes_disk.dat
		awk 'BEGIN { print int((${tBamSize}+${tBaiSize}+${nBamSize}+${nBaiSize}+${dbSNPVCFSize}+${cosmicVCFSize}+${normalPanelSize}+${fastaSize}+${fastaDictSize}+3000000000)/1000000000+1) }' > $SIZE_FILE

		#Create list of indices for the scatter job
		seq 0 $((${nWay}-1)) > indices.dat

		#Run the prepare task that splits the .interval_list file into subfiles
		MUTECT_INTERVALS="${if defined(sampleTargetsIntervalList) then sampleTargetsIntervalList else workspaceTargetsIntervalList}"
		java -Xmx2g -jar /usr/local/bin/GatkScatterGatherPrepare.jar . ${nWay} --intervals $MUTECT_INTERVALS --reference_sequence ${refFasta}

		>>>

	output  {
		Array[File] interval_files=glob("gatk-scatter*")
		Int mutectDiskGB=read_int("split_base_sizes_disk.dat")
		Array[Int] scatterIndices=read_lines("indices.dat")
		}

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 100 HDD"
        continueOnReturnCode: [0, 1]
		}
	}
task Mutect1_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible
		Int diskGB
		#TASK INPUT PARAMS
		File mutectIntervals
		File tumorBam
		#File tumorBamIdx
		File normal_Bam
		File normal_BamIdx
		File refFasta
		File refFastaDict
		File refFastaIdx
		Float fracContam
		Float contamFloor
		File dbSNPVCF
    File dbSNPVcfIndex
		File cosmicVCF
		Int downsampleToCoverage
		File? readGroupBlackList
		File? normalPanel
        String normalBam
        String normalBamIdx
		String pairName
		String caseName
		String bam_current_dir
		String bai_current_dir
		String? ctrlName

		String? optional1
		String? optional2
		String? optional3
		String? optional4
		String? optional5
	command <<<
	#increase verbosity
	set -x

	#convert tumor bam file to bai file using samtools

	cp ${tumorBam} .
	samtools index -b ${tumorBam}
	cp ${tumorBam}.bai .
	
    # Additions
	samtools view -H ${normal_Bam} | cat  
	
    
    ln -s "${normal_Bam}" "${normalBam}.cram"
	ln -s "${normal_BamIdx}" "${normalBamIdx}.cram.crai"
	
    
    
    #compute/apply contamination floor for effective contamination
	EFFECTIVE_CONTAMINATION=`/usr/local/bin/python -c 'import sys;print sys.argv[1] if(float(sys.argv[1])>=float(sys.argv[2])) else sys.argv[2]'  ${fracContam} ${contamFloor}`
	echo "EFFECTIVE_CONTAMINATION computed to be $EFFECTIVE_CONTAMINATION"

	#mutect 1
	java -jar -Xmx9g /usr/local/bin/muTect-1.1.6.jar --analysis_type MuTect \
	-L ${mutectIntervals} \
    -U ALLOW_N_CIGAR_READS \
	${if defined(ctrlName) then "${ctrlName} " else " "} \
	--tumor_sample_name ${caseName} -I:tumor ${tumorBam} \
	--reference_sequence ${refFasta} \
	--fraction_contamination $EFFECTIVE_CONTAMINATION --dbsnp ${dbSNPVCF} \
	--cosmic ${cosmicVCF} ${"--read_group_black_list" + readGroupBlackList}\
	--out ${pairName}.call_stats.txt --coverage_file ${pairName}.coverage.wig.txt \
	--power_file ${pairName}.power.wig.txt  --downsample_to_coverage ${downsampleToCoverage} \
	--enable_extended_output ${optional1} ${optional2} ${optional3} ${optional4} ${optional5} \
    ${if defined(normalPanel) then "${normalPanel} " else " "} \
	>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		memory: "12 GB"
		disks: "local-disk 100 HDD"
		}


	output {
		File tumorBamIdx="${tumorBam}.bai"
		File mutect1_cs="${pairName}.call_stats.txt"
		File mutect1_pw="${pairName}.power.wig.txt"
		File mutect1_cw="${pairName}.coverage.wig.txt"
		}

	}
task Gather_Task {
		#RUNTIME INPUT PARAMS
		Int preemptible

		#TASK INPUT PARAMS
		Array[File] mutect1_cs
		Array[File] mutect1_pw
		Array[File] mutect1_cw
		String pairName
	command <<<
		#increase verbosity
		set -x

		#mutect1 call_stats merging
		MUTECT1_CS="${pairName}.call_stats.txt"
		head --lines=2 ${mutect1_cs[0]} > $MUTECT1_CS
		cat ${sep =' ' mutect1_cs} | grep -Pv '#'|grep -Pv '^contig' >> $MUTECT1_CS

		MUTECT1_PW="${pairName}.power.wig.txt"
		MUTECT1_PW_GZ="${pairName}.power.wig.txt.gz"
		head --lines=1 ${mutect1_pw[0]} > $MUTECT1_PW
		cat ${sep =' ' mutect1_pw} | grep -Pv '^track' >> $MUTECT1_PW
		tar -zcvf $MUTECT1_PW_GZ $MUTECT1_PW

		MUTECT1_CW="${pairName}.coverage.wig.txt"
		MUTECT1_CW_GZ="${pairName}.coverage.wig.txt.gz"
		head --lines=1 ${mutect1_cw[0]} > $MUTECT1_CW
		cat ${sep =' ' mutect1_cw} | grep -Pv '^track' >> $MUTECT1_CW
		tar -zcvf $MUTECT1_CW_GZ $MUTECT1_CW
		>>>


	output {
		File Mutect_call_stats="${pairName}.call_stats.txt"
		File Mutect_power_wig="${pairName}.power.wig.txt.gz"
		File Mutect_coverage_wig="${pairName}.coverage.wig.txt.gz"
		}


	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 100 HDD"
		}
	}

#This task corresponds to CallstatsToMaflite_V14 task under the "Call Somatic Mutations For Capture" workflow.
task CallstatsToMaflite_14 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File callstatsFile
		Int genomeBuild
		String mode
		String pairName
		String extraColumns
		String outputFile = "${pairName}.maf"
	command <<<
		#increase verbosity
		set -x
		echo "Calling call_stats_to_maflite.pl on input file ${callstatsFile}. Dumping to ${outputFile}"

		perl /usr/local/bin/call_stats_to_maflite.pl ${callstatsFile} ${genomeBuild} ${mode} ${outputFile} "${extraColumns}"

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 100 HDD"
		}

	output {
		File maflite_capture="${outputFile}"
	}
}

task SummarizeWigFile_1 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File coverageFile
		String pairName
		String outputFile = "${pairName}.somatic_coverage_summary.txt"
	command <<<
		#increase verbosity
		set -x
		echo "Calling summarizeWigFile.pl on input file ${coverageFile}. Dumping to ${outputFile}"


		tar -zxvf ${coverageFile}
		perl /usr/local/bin/summarizeWigFile.pl "${pairName}.coverage.wig.txt" ${outputFile}

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 100 HDD"
		}

	output {
		File somatic_mutations_covered_bases_capture="${outputFile}"
	}
}

task ApplyMafValidation_23 {
		#RUNTIME INPUT PARAMS
		Int preemptible

		File mafliteFile
		String sampleName
		String matchMode
		File valdb
		String valDir = sub(valdb, "\\.tar.gz$", "")

		String outputFile="${sampleName}.maf.annotated"
	command <<<
		#increase verbosity
		set -x
		echo "Calling ApplyMafValidation.java on input file ${mafliteFile}. Dumping to ${outputFile}"

		VAL_BASENAME=`basename ${valdb}`
		VAL_DIR=`echo $VAL_BASENAME | cut -d'.' -f1`
		tar -zxvf ${valdb}

		java -Xmx1g -jar /usr/local/bin/ApplyMAFValidation.jar M=${mafliteFile} OUTPUT_MAF=${outputFile} MATCH_MODE=${matchMode} V=$VAL_DIR

		>>>

	runtime {
		preemptible: "${preemptible}"
		docker: "bbecerra/rna_mutect:scripts_added"
		disks: "local-disk 100 HDD"
		}

	output {
		File snp_validated_maflite_capture="${outputFile}"
	}
}

workflow CallingGroup_Workflow {
		#RUNTIME INPUT PARAMS
		Int preemptible

		#WORKFLOW INPUT PARMS
		#General inputs for the whole workflow
		String pairName
		String tumorName
		File tumorBamRNA
		#File tumorBamRNAIdx
		String? normalName
		File normalBamDNA
		File normalBamDNAIdx
		File refFasta
		File refFastaDict
		File refFastaIdx

		#Input for BamRealignment
		File originalRefFasta
		File originalRefFastaIndex

		File refBwt
		File refSa
		File refAmb
		File refAnn
		File refPac
		File dbSNPVcfIndex

		Array[File] knownIndelsSitesVCFs
		Array[File] knownIndelsSitesIndices

		#Inputs for CallSomaticMutations_131_Prepare_Task
		Int scatterNWay
		File? sampleTargetsIntervalList
		File workspaceTargetsIntervalList

        #Inputs for MasterFilter_4
        File? snpFilteringPoNFile

		#Inputs for Mutect1_Task
		Float fracContam
		File dbSNPVCF
		File cosmicVCF
		Int downsampleToCoverage
		File? readGroupBlackList
		File? normalPanelVCF
		String? optional1
		String? optional2
		String? optional3
		String? optional4
		String? optional5

		#Inputs for CallstatsToMaflite
		Int genomeBuild_int
		String mode
		String extraColumns

		# Inputs for ApplyMafValidation
		File valdb
		String matchMode

		#Inputs for Oncotator
		String genomeBuild_string
		File oncoDBTarBall
		File defaultConfig
		String inputFormat
		String outputFormat
		String txMode
		String additionalOptions
		File? canonicalTXFile
		String? logFilename


		#Inputs for HiSat2Realign_3
		File HISAT_index
		Int? min_intronlen
		Int? max_intronlen
		File? known_splicesite_infile
		File? novel_splicesite_infile
		String? rna_strandness
		Int? k
		Int? minins
		Int? maxins
		String? mate_orientation
		String? rg_id
		String? rg
		Int? num_threads

		#Inputs for CallSomaticMutations_131_Prepare_Task after hisat realignment
		Int hisat_scatterNWay
		File? hisat_sampleTargetsIntervalList

		#Inputs for Mutect1_Task after hisat realignment
		Int hisat_downsampleToCoverage
		File? hisat_readGroupBlackList
		String? hisat_optional1
		String? hisat_optional2
		String? hisat_optional3
		String? hisat_optional4
		String? hisat_optional5

		#Inputs for FilterRNAMutations
		File DarnedMat
		File ExacMat
		File RadarMat
        File? ponGTEx
		Int minAlt_filter_generic
        Int? PoNThr
	    File cytoBand
		#Inputs for FilterRNAMutations_duplicateReads
		Int minAlt_filter_duplicates

		String samtools_docker = "jweinstk/samtools:latest"

        Int m2BufferGB = 15

	# PREPARE FOR SCATTER
	call CallSomaticMutations_131_Prepare_Task as CallSomaticMutations_131_Prepare_Task_STAR{
		input:
			refFasta=refFasta,
			refFastaDict=refFastaDict,
			refFastaIdx=refFastaIdx,
			sampleTargetsIntervalList=sampleTargetsIntervalList,
			workspaceTargetsIntervalList=workspaceTargetsIntervalList,
			nWay=scatterNWay,
			preemptible=preemptible,
			tBamSize=size(tumorBamRNA),
			#tBaiSize=size(tumorBamRNAIdx),
			#nBamSize=size(normalBamDNA),
			#nBaiSize=size(normalBamDNAIdx),
			dbSNPVCFSize=size(dbSNPVCF),
			cosmicVCFSize=size(cosmicVCF),
			#normalPanelSize=size(normalPanelVCF),
			fastaSize=size(refFasta),
			fastaDictSize=size(refFastaDict)

	}

	#SCATTER AND ANALYZE
	scatter (idx in CallSomaticMutations_131_Prepare_Task_STAR.scatterIndices) {

			call Mutect1_Task as Mutect1_Task_STAR {
				input:
					contamFloor=0.01,
					tumorBam=tumorBamRNA,
					normal_Bam=normalBamDNA,
					#tumorBamIdx=tumorBamRNAIdx,
					pairName=pairName,
					caseName=tumorName,
					#ctrlName=normalName,
					normal_BamIdx=normalBamDNAIdx,
					mutectIntervals=CallSomaticMutations_131_Prepare_Task_STAR.interval_files[idx],
					refFasta=refFasta,
					refFastaDict=refFastaDict,
					refFastaIdx=refFastaIdx,
					fracContam=fracContam,
					dbSNPVCF=dbSNPVCF,
                    dbSNPVcfIndex=dbSNPVcfIndex,
					cosmicVCF=cosmicVCF,
					downsampleToCoverage=downsampleToCoverage,
					readGroupBlackList=readGroupBlackList,
					#normalPanel=normalPanelVCF,
					preemptible=preemptible,
					diskGB=CallSomaticMutations_131_Prepare_Task_STAR.mutectDiskGB + m2BufferGB,
					optional1=optional1,
					optional2=optional2,
					optional3=optional3,
					optional4=optional4,
					optional5=optional5

				}


			}


	call Gather_Task as Gather_Task_STAR {
		input:
			pairName=pairName,
			mutect1_cs=Mutect1_Task_STAR.mutect1_cs,
			mutect1_pw=Mutect1_Task_STAR.mutect1_pw,
			mutect1_cw=Mutect1_Task_STAR.mutect1_cw,
			preemptible=preemptible,

	}


	call CallstatsToMaflite_14 {
		input:
			preemptible=preemptible,
			callstatsFile=Gather_Task_STAR.Mutect_call_stats,
			genomeBuild=genomeBuild_int,
			mode=mode,
			pairName=pairName,
			extraColumns=extraColumns
	}

	call SummarizeWigFile_1 {
		input:
			preemptible=preemptible,
			coverageFile=Gather_Task_STAR.Mutect_coverage_wig,
			pairName=pairName
	}

	call ApplyMafValidation_23 {
		input:
			preemptible=preemptible,
			mafliteFile=CallstatsToMaflite_14.maflite_capture,
			sampleName=tumorName,
			valdb=valdb,
			matchMode=matchMode
	}

	output {
		File Mutect_call_stats_STAR = Gather_Task_STAR.Mutect_call_stats
		File Mutect_power_wig_STAR = Gather_Task_STAR.Mutect_power_wig
		File Mutect_coverage_wig_STAR = Gather_Task_STAR.Mutect_coverage_wig
		File maflite_capture_STAR = CallstatsToMaflite_14.maflite_capture
		File somatic_mutations_covered_bases_capture_STAR = SummarizeWigFile_1.somatic_mutations_covered_bases_capture
		File snp_validated_maflite_capture_STAR = ApplyMafValidation_23.snp_validated_maflite_capture
        }
}
