#python libraries
import subprocess
import sys
import re

#third party libraries
import abifpy

# biopython goodies
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO

class SeqFixer():
    def __init__(self,filename,cdrRecords,reportFile):
        self.filename     = filename
        self.cdrRecords   = cdrRecords
        self.reportFile   = reportFile
        self.reportHandle = open(reportFile,'a')


        full_file = self.filename.split('/')[-1]
        root_name = full_file.split('.')[0]

        self.reportHandle.write('In sequence file '+ root_name + ' the fixes made are:\n\n')
        self.cdrModifications = False
 
    def callBases(self):
        record = self.rawRecord
        
        badCallSeq   = [ind for ind, base  in enumerate(record.seq) if base=='N']    
        badCallTrace = [peak for base, peak in zip(record.seq,record.data['tracepeaks']) if base=='N']

        gTrace = record.tags['DATA9'].tag_data
        aTrace = record.tags['DATA10'].tag_data
        tTrace = record.tags['DATA11'].tag_data
        cTrace = record.tags['DATA12'].tag_data

        traceSums = [
            [sum(gTrace[(i-2):(i+2)]),
             sum(aTrace[(i-2):(i+2)]),
             sum(tTrace[(i-2):(i+2)]),
             sum(cTrace[(i-2):(i+2)])]
            for i in badCallTrace]

        revisedCalls = [sums.index(max(sums))   for sums in traceSums]
        revisedSeq = list(record.seq)

        baseOrder = ['G','A','T','C']
        for pos,baseInd in zip(badCallSeq,revisedCalls):
            revisedSeq[pos] = baseOrder[baseInd]
            self.reportHandle.write('Base at position {0!s} was recalled as {1}\n'.format(pos,baseOrder[baseInd]))
        self.reportHandle.write('\n')
        
        # enter report info on revised calls
        
        self.seqWithBasesCalled = SeqRecord(Seq(''.join(revisedSeq),DNAAlphabet),id='BlankID')

    def cdrAlign(self,records):

        # initialize muscle commandline call
        muscle_cline = MuscleCommandline()
        child = subprocess.Popen(str(muscle_cline),
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=(sys.platform!="win32"))
        SeqIO.write(records,child.stdin,"fasta")
        child.stdin.close() # must close handle to start process

        self.alignedRecords = AlignIO.read(child.stdout,"fasta")
        # print self.alignedRecords[0].seq
        # print self.alignedRecords[1].seq

    def fixGeneral(self):
            
        # main logic here
        bases = re.compile('G|A|T|C')

        # find boundaries

        match_obj = bases.search(self.alignedRecords[0].seq.tostring())
        start_index = match_obj.start()
        cdr_len = len(bases.findall(self.alignedRecords[0].seq.tostring()))

        # copy in original cdr

        localFixedRecord = ''.join([
            self.alignedRecords[1].seq.tostring()[0:start_index],
            self.alignedRecords[0].seq.tostring()[(start_index):(start_index+cdr_len)],
            self.alignedRecords[1].seq.tostring()[(start_index+cdr_len):]
            ])

        # remove any '-' from sequence. Possible if alignment created '-' in CDR
        localFixedRecord = SeqRecord(Seq(localFixedRecord.replace('-',''),DNAAlphabet),id="fixedSequence")
        if not (localFixedRecord.seq.tostring() == self.fixedRecord.seq.tostring()):
            self.reportHandle.write('Using alignment between CDR and sequence:\n')
            self.reportHandle.write(self.alignedRecords[0].seq.tostring())
            self.reportHandle.write('\n')
            self.reportHandle.write(self.alignedRecords[1].seq.tostring())
            self.reportHandle.write('\n')
            self.reportHandle.write('The sequence was corrected to be:\n')
            self.reportHandle.write(localFixedRecord.seq.tostring())
            self.reportHandle.write('\n')

            self.cdrModifications = True
            
        self.fixedRecord = localFixedRecord

        # check if fixedRecord is same as input, if not, record info about fix
        
    # top function call
    def fixSequence(self):

        # define in memory SeqRecords
        self.rawRecord = abifpy.Trace(self.filename)

        # require unambiguous alphabet
        self.callBases()
        self.fixedRecord = self.seqWithBasesCalled

         # multiple copies of the same CDR would be a problem
        for cdrRec in self.cdrRecords:

            records = (cdrRec,self.fixedRecord)

            self.cdrAlign(records)
            self.fixGeneral()
        
        return self.fixedRecord


