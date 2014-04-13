
#!/usr/bin/python
import sys, os, re, getopt, glob, numpy as np
import timeit
import itertools

start = timeit.default_timer()

Usage = """
Usage:  ./grab_ALL_READS_from_mothur_list_file.py -c 0.01 -l mothur_list_file.an.list -e Fasta_enr_otu.csv

REQUIRED ARGUMENTS:

                -c      the sequence similarity cut-off used in creating the OTU definition

                -l      the list file produced in the final stages of the mothur pipeline, which group all reads into OTUs

                -e      list of all enrOTU read names, one per line.
		
		-f 	fasta file containing all reads (i.e. unprocessed by mothur)

OPTIONAL:
                -b <V>  Enter DEBUG MODE - saves debugging information to "error.log" 
                        (Must proivde "V")

		-R	Include taxonomic affiliations by using output from R script entitled:  "Block05 - Converting List of Read Names into Taxonomic Descriptions + Read Names"

		-m	Process the reads with pre-set mothur commands (change them within the script if ye desire)
			(Needs .qual and .oligos files + you must manually change the names in the body of this script)
UTILITY:
	This script performs two key functions in SIP-DNA:
		i) It pulls ALL fasta sequences from within each enrOTUs (at a user specified definitions of OTU).
		ii) It pulls sequence data for all enrOTU reps
"""

if len(sys.argv)<2:
        print Usage
        sys.exit()


# Store input and output file names
CUT_OFF=''
LIST_FILE=''
ENR_OTUs=''
FASTA=''
TAXA=''
MOTHUR=''
Debug=''


# Read command line args
myopts, args = getopt.getopt(sys.argv[1:],"l:e:c:b:f:R:m:")

###############################
# o == option
# a == argument passed to the o
###############################
for o, a in myopts:
    if o == '-c':
        CUT_OFF= a
    if o == '-l':
        LIST_FILE= a
    if o == '-e':
        ENR_OTUs= a
    if o == '-f':
        FASTA= a
    if o == '-R':
        TAXA= a
    if o == '-m':
        MOTHUR= a
    if o == '-b':
        Debug= a

if Debug:
        print "You are in Debugging Mode"
        error = open("error.log", "w")

#Extract OTU List According to Desired Cut-off
END = "Haven't Even Begun"
while END != "YES":
	for line in open(LIST_FILE, "r"):
        	        line = line.split()

			#Select line which contains cut-off
			if re.search(CUT_OFF, line[0]):
				LINE = line
				END = "YES"
if Debug:
        error.write("Your List File Says there Are: "+str(LINE[1])+" OTU Groups. We're Processing: "+str(np.size(LINE)-2)+". Make sure these numbers match.\n")
	error.write("Did We Find your Cut-off? "+ END+"\n")


#Import Taxa Details if Specified
TAXA_DICT={}

if TAXA: 
	with open(TAXA, "r") as infile:
	        for line in infile:
        	        line = line.split()
			if re.search(">",line[0]):
	                        READ_NAME = re.sub(">", "", line[0])
        	                TAXA_DICT[READ_NAME] = ' '.join(line[1:np.size(line)])


#Make Two Dictionaries; one referencing each read to OTU group and the other vice versa
OTU_DICT={}	#READNAME: OTU_Number
GROUP_DICT={}	#OTU_Number: READ_NAME1, READ_NAME2 ... READ_NAMEn

for ante in range(np.size(LINE)):
	if ante > 1:
		GROUP = LINE[ante]
		GROUP = GROUP.split(",")

		if np.size(GROUP) > 1:
			for OTU in GROUP:
				OTU_DICT[OTU] = ante

				if not GROUP_DICT.has_key(ante):
					GROUP_DICT[ante] = []

                        	GROUP_DICT[ante].append(OTU)

		else:
			OTU_DICT[OTU] = ante
			GROUP_DICT[ante] = GROUP

#Make Dictionary of FASTA Reads
FASTA_DICT = {}

with open(FASTA, "r") as infile:
	for line in infile:
		line = line.split()
		if re.search(">",line[0]):
			header = re.sub(">", "", line[0])
			FASTA_DICT[header] = next(infile)


#Create Output File
output = open("./Output_enrOTU_only.fasta", "w")			#FASTA File of enrOTU Only
output2 = open("./Output_All_READS_from_enrOTU.fasta", "w")		#LIST OF original enrOTU + Read name
output3 = open("./Output_LIST_of_All_READS_from_enrOTU.list", "w")	#FASTA File of All Reads

#Grab enrOTUs
for line in open(ENR_OTUs, "r"):
		enrOTU = line.strip("\n\r")		
		
		if Debug:
			error.write(str(len(enrOTU))+"\n")

		if len(enrOTU) == 0:
			pass

		else:
			if TAXA:
				#Grab Taxa Information
				Long_Header_Info = TAXA_DICT[enrOTU]

				if Debug:
					error.write(Long_Header_Info)

				#Print enrOTU only FASTA File				
				output.write(">"+enrOTU+" "+Long_Header_Info+"\n"+FASTA_DICT[enrOTU])

				#Find OTU GROUP
				Group_ID = OTU_DICT[enrOTU]

				if Debug:
					error.write(str(Group_ID)+"\n")
	
				#Grab all members of Group
				Total_Group = GROUP_DICT[Group_ID]

				#Print ALL READs FASTA File
				for printme in Total_Group:
					output2.write(">"+printme+" "+Long_Header_Info+"\n"+FASTA_DICT[printme])
					output3.write(printme+"_enrOTUrepseq_"+Total_Group[0]+"\n")

			else:
				#Print enrOTU only FASTA File				
				output.write(">"+enrOTU+"\n"+FASTA_DICT[enrOTU])

				#Find OTU GROUP
				Group_ID = OTU_DICT[enrOTU]

				if Debug:
					error.write(str(Group_ID)+"\n")
	
				#Grab all members of Group
				Total_Group = GROUP_DICT[Group_ID]

				#Print ALL READs FASTA File and List
				for printme in Total_Group:
					output2.write(">"+printme+"\n"+FASTA_DICT[printme])
					output3.write(printme+"_enrOTUrepseq_"+Total_Group[0]+"\n")

output.close()
output2.close()
output3.close()

if MOTHUR:
	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.fasta; do mothur \"# trim.seqs(fasta=$n, oligos=IIKFCBR01.oligos, qfile=IIKFCBR01.qual, maxambig=0, maxhomop=8, flip=T, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, processors=10)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.fasta; do mothur \"# unique.seqs(fasta=$n)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.fasta; do mothur \"# align.seqs(fasta=$n, reference=~/Phylogenetic_Gene_Databases/silva.bacteria/silva.bacteria.fasta, flip=T, processors=12)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.align; do mothur \"# screen.seqs(fasta=$n, name=Output_All_READS_from_enrOTU.trim.names, start=1044, minlength=250, processors=12)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.good.align; do mothur \"# filter.seqs(fasta=$n, vertical=T, processors=12)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.good.filter.fasta; do mothur \"# unique.seqs(fasta=$n, name=Output_All_READS_from_enrOTU.trim.good.names)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.good.filter.unique.fasta; do mothur \"# chimera.uchime(fasta=$n, name=Output_All_READS_from_enrOTU.trim.unique.good.filter.names, processors=10)\"; done"
	]))

	os.system(' '.join([
		"for n in Output_All_READS_from_enrOTU.trim.unique.good.filter.unique.uchime.accnos; do mothur \"# remove.seqs(accnos=$n, fasta=Output_All_READS_from_enrOTU.trim.unique.good.filter.unique.fasta, name=Output_All_READS_from_enrOTU.trim.unique.good.filter.names)\"; done"
	]))

	os.system(' '.join([
		"mv Output_All_READS_from_enrOTU.trim.unique.good.filter.pick.names enrOTU_final.names"
	]))

	os.system(' '.join([
		"mv Output_All_READS_from_enrOTU.trim.unique.good.filter.unique.pick.fasta enrOTU_final.fasta"
	]))

stop = timeit.default_timer()

print stop - start

