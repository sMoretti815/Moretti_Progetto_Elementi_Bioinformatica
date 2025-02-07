#!/usr/bin/env python3

import argparse
import pysam
from pysam import AlignmentFile
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--fileGtf", default = "./annotation_one_tr_chr21.gtf", type = argparse.FileType('r'), help="Nome del file gtf")
parser.add_argument("-b", "--fileBam", default= "./sample-chr21.bam", type = argparse.FileType('rb'), help="Nome del file Bam")
parser.add_argument("-o", "--output", default = "spliced_reads", type = str, help="Nome del file di output")
args = parser.parse_args()

input_bam = args.fileBam
pysam.index(input_bam.name)
bam_file = AlignmentFile(input_bam, "rb")

input_gtf = args.fileGtf
df = pd.read_csv(input_gtf, sep="\t", header=None)
replace_dict = {0 : 'reference', 1 : 'source', 2 : 'feature', 3 : 'start',
                4 : 'end', 5 : 'score', 6 : 'strand', 7 : 'frame', 8 : 'attributes'}
df.rename(columns = replace_dict, inplace = True)
df["transcript_id"] = df["attributes"].apply(lambda x : re.search(r"transcript_id\s"+ "(.+?);", x).group(1))
df["gene_id"] = df["attributes"].apply(lambda x : re.search(r"gene_id\s" + "(.+?);", x).group(1))

if bam_file.count() == 0:
    print("il file bam non contiene reads")
    exit(0)
    
if (bam_file.nreferences > 1):
    print("il file bam contiene più di una reference")
elif (bam_file.mapped == 0):
    print("il file non contiene reads mappati")
elif len(df.reference.unique()) > 1:
    print("il file gtf contiene più di una reference")
else:
    if str(df.reference.unique()[0]) != str(bam_file.references[0]):
        print("la reference del file gtf potrebbe non corrispondere con quella del file bam")
    
    start_mapped_position = float("inf")
    end_mapped_position = float("-inf")
    for read in bam_file.fetch():
        if read.is_mapped:
            start_mapped_position = min(start_mapped_position, read.reference_start)
            end_mapped_position = max(end_mapped_position, read.reference_end)
    start_mapped_position += 1
    end_mapped_position += 1
    print("tutti i reads nel file bam sono mappati tra le posizioni " + str(start_mapped_position)
           + " e " + str(end_mapped_position) + " della reference")

    df_grouped = df[["start", "end", "transcript_id", "gene_id"]].groupby(["transcript_id","gene_id"]).agg({"start": "min", "end": "max"})
    df_grouped.reset_index(inplace = True)

    max_position_gtf = max(df_grouped.end) 
    min_position_gtf = min(df_grouped.start)

    if start_mapped_position < min_position_gtf or end_mapped_position > max_position_gtf:
        print("i reads sono stati sequenziati fuori dal range coperto dal file gtf")
    else:
        df_grouped = df_grouped[(df_grouped.start <= start_mapped_position) & (df_grouped.end >= end_mapped_position)]

        if df_grouped.empty:
            print("i reads sono stati sequenziati da trascritti diversi")
        else:
            transcript_possibili = list(zip(df_grouped.transcript_id, df_grouped.gene_id))
            
            if len(transcript_possibili) > 1: 
                df = df[(df.end >= start_mapped_position) & (df.start <= end_mapped_position) & (df.feature == "exon")]
                df.reset_index(inplace = True)
                introns = sorted(bam_file.find_introns(bam_file.fetch()))
                for (transcript_id, _) in list(transcript_possibili):
                    if len(df[transcript_id == df["transcript_id"]]) != len(introns)+1:
                        transcript_possibili.remove((transcript_id, _))
                    else:
                        for (i, (start_intron, end_intron)) in enumerate(introns, start=df[df["transcript_id"] == transcript_id].index[0]):
                            if (df.loc[i].end != start_intron or df.loc[i+1].start != end_intron+1):
                                transcript_possibili.remove((transcript_id, _))
                                break
        
            print("i reads sono stati sequenziati dal trascritto " + transcript_possibili[0][0] + " del gene " + transcript_possibili[0][1])

with open(args.output + ".fastq", "w") as output_file:
    for read in bam_file.fetch():
        if "N" in read.cigarstring:
            qualities = read.query_qualities if read.query_qualities != None else [0]*len(read.query_sequence)
            record = SeqRecord(Seq(read.query_sequence), id=read.query_name,
                               description = "base " + str(read.next_reference_start) + " on " + read.reference_name,
                               letter_annotations={"phred_quality": qualities})
            output_file.write(format(record,"fastq"))