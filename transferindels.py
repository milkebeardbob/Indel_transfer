#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2025 Julius Guse & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script transfers the gaps from a set of alignments (FASTA format) to corresponding (same names, same sequence accessions) unaligned FASTA files. It is meant to transfer the gaps to reduced amino acid alphabet datasets.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that works on all files in the working directory with a given extension.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)

import sys
import os

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio" : "https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from Bio import SeqIO
from Bio.Seq import Seq

print('#Script: transferindels.py')
print('#Version: v2025_1')
print('#Usage: python transferindels.py <unaligned> <unaligned_ext> <alignments> <alignments_ext> <output_ext>')
print('#<unaligned> must be the directory containing the unaligned sequence FASTA files with <unaligned_ext>. (trailing slash optional) (required)')
print("#<unaligned_ext> must be the filename extension of the unaligned sequence FASTA files. (required)")
print('#<alignments> must be the directory containing the aligned sequence FASTA files with <alignments_ext>. (trailing slash optional) (required)')
print("#<alignments_ext> must be the filename extension of the aligned sequence FASTA files. (required)")
print('#<output_ext> must be the extension of the output FASTA files with transferred indels. (leading dot optional) (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

#Check if the correct number of arguments is given
if len(sys.argv) == 6:
    print ('Five arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

#We don't assign variables to the directory paths now, because they will be converted to an abspath later.
uaext = sys.argv[2]
alext = sys.argv[4]
outext = sys.argv[5]

#Check if the extensions start with a dot, otherwise add them.
if uaext.startswith('.') == False: #This is not absolutely necessary, since os.path.splitext will detect the extension anyway. It's more of a precaution against double dots.
    uaext = str('.' + uaext)
if alext.startswith('.') == False:
    alext = str('.' + alext)
if outext.startswith('.') == False:
    outext = str('.' + outext)

#Checkpoint for datasets directory existence and trailing slash. Convert to abspath to make sure there are no issues when called through doggo_zoomies.
if os.path.exists(sys.argv[1]) == True:
    print ('Unaligned dataset directory found. Proceeding.')
    uadir = os.path.abspath(sys.argv[1])
    uadir = os.path.join(uadir, '')
else:
    print ('Unaligned dataset directory not found. Exiting.')
    sys.exit(1)

if os.path.exists(sys.argv[3]) == True:
    print ('Aligned dataset directory found. Proceeding.')
    aldir = os.path.abspath(sys.argv[3])
    aldir = os.path.join(aldir, '')
else:
    print ('Aligned dataset directory not found. Exiting.')
    sys.exit(1)

#Check if files with a given extension exist in both the unaligned and aligned datasets directories and create a list of each.
uafilenames = []
for uafname in os.listdir(uadir):
    if uafname.endswith(uaext):
        uafname = os.path.join(uadir, uafname)
        uafilenames.append(uafname)
if len(uafilenames) > 0:
    print('File(s) with the unaligned dataset extension found in the unaligned dataset directory. Proceeding.')
else:
    print('No files with the unaligned dataset extension found in the unaligned dataset directory. Exiting.')
    sys.exit(1)

alfilenames = []
for alfname in os.listdir(aldir):
    if alfname.endswith(alext):
        alfname = os.path.join(aldir, alfname)
        alfilenames.append(alfname)
if len(alfilenames) > 0:
    print('File(s) with the aligned dataset extension found in the aligned dataset directory. Proceeding.')
else:
    print('No files with the aligned dataset extension found in the aligned dataset directory. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print('Removing files with names identical to the output.')
removal = ('rm -r *' + outext + ' 2> /dev/null')
os.system(removal)

#Iterate over each file in the list for unaligned datasets.
for iteruafname in uafilenames:
    corralfname = str(aldir + os.path.basename(iteruafname).split(os.extsep, 1)[0] + alext) #Find the corresponding aligned dataset, with which they share a file stem.
    if not os.path.exists(corralfname):
        print ('WARNING: Couldn\'t find aligned dataset corresponding to ' + iteruafname)
        #Decided not to exit with error here just in case there are different datasets being matched. TODO: Exit instead. Checkpoint should be placed when checking for files present in the directories. Check should be for existence of aligned corresponding to aligned, since the gaps will be transferred there. Vice versa is pointless.
        continue

    # Read unaligned dataset sequences
    uarecords = SeqIO.to_dict(SeqIO.parse(iteruafname, "fasta"))

    # Read corresponding aligned dataset sequences
    alrecords = SeqIO.to_dict(SeqIO.parse(corralfname, "fasta"))

    # Create new empty list of SeqRecords with indels transferred.
    updated_records = []

    # Ensure that the accessions from each dataset are present in the other.
    for acc1 in uarecords:
        if acc1 not in alrecords:
            print('Accession ' + acc1 + ' from unaligned dataset ' + iteruafname + ' not found in aligned dataset ' + corralfname + '. Exiting.')
            sys.exit(1)

    for accession in alrecords:
        if accession not in uarecords:
            print('Accession ' + accession + ' from aligned dataset ' + corralfname + ' not found in unaligned dataset ' + iteruafname + '. Exiting.')
            sys.exit(1)

        aligned_al_seq = alrecords[accession].seq #Care that the seq gets assigned to a string here, see comments below.
        original_ua_seq = uarecords[accession].seq

        #Check that unaligned and aligned sequences have the same number of residues.
        if aligned_al_seq.count("-") + len(original_ua_seq) != len(aligned_al_seq):
            print('Unaligned and aligned sequences for ' + accession + ' in dataset ' + iteruafname + ' do not have the same number of residues. Exiting.')
            sys.exit(1)

        new_ua_seq = [] #Start new sequence as list and index.
        ua_index = 0

        #If you see a gap in aligned, add gap. If you see a residue, add the corresponding residue number (based on index) of unaligned and increase index by 1. We assume that the sequences are the same and individual residues have not been changed somehow. TODO: Add a check?
        for res in aligned_al_seq:
            if res == "-":
                new_ua_seq.append("-")
            else:
                new_ua_seq.append(original_ua_seq[ua_index])
                ua_index += 1

        updated_record = alrecords[accession]
        updated_record.seq = Seq("".join(new_ua_seq)) #This syntax because as of Python 3.12: "BiopythonDeprecationWarning: Using a string as the sequence is deprecated and will raise a TypeError in future. It has been converted to a Seq object."
        updated_records.append(updated_record)

    # Write output to new file with outext extension.
    outalname = str(os.path.basename(iteruafname).split(os.extsep, 1)[0] + outext)
    with open(outalname, "w") as out_f:
        SeqIO.write(updated_records, out_f, "fasta")

print('All done!')
