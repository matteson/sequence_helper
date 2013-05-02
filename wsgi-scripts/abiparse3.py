import abifpy
import glob

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from SeqFixerClass import SeqFixer

def abiparse(directory, glob_pattern,cdr_file):

    records_summary = []

    cdr_mods_list_handle = open(directory + 'CDR_modifications.txt','a')

    for filename in glob.glob(directory + glob_pattern):
        cdr_input_handle = open(cdr_file,"rU")

        full_file = filename.split('/')[-1]
        print full_file
        root_name = full_file.split('.')[0]
        # open and parse the abi file, shoul do some error handling here

        # v1f1 = SeqRecord(Seq('gatatccagatgacccagtccccctcgagcctgtccgcctctgtgggcgatagggtcaccatcacctgccgtgccagtcag',DNAAlphabet),id="testcdr")

        test = SeqFixer(filename,list(SeqIO.parse(cdr_input_handle,"fasta")),directory + root_name + '_report.txt')

        out = test.fixSequence()
        test.reportHandle.close()

        # write fasta files
        fasta_out = open( directory + root_name + '.fasta','w')
        # test.fixedRecord.id = root_name
        SeqIO.write(test.fixedRecord, fasta_out,'fasta')
        fasta_out.close()

        records_summary.append(test.fixedRecord)

        if test.cdrModifications:
            cdr_mods_list_handle.write('CDR modifications made for sequence: ' + root_name + '\n')
            print 'this should be in the cdr list: ' + root_name    

        cdr_input_handle.close()

    fasta_summary_out = open(directory + 'fasta_summary.fasta', 'w')
    SeqIO.write(records_summary,fasta_summary_out,'fasta')
    fasta_summary_out.close()

    cdr_mods_list_handle.close()
