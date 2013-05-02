import abifpy
import glob
from Bio import SeqIO

def abiparse(directory, glob_pattern):

    for filename in glob.glob(directory + glob_pattern):
        full_file = filename.split('/')[-1]
        print full_file
        root_name = full_file.split('.')[0]
        # open and parse the abi file, shoul do some error handling here
        test = abifpy.Trace(filename)

        handle_out = open( directory + root_name + '.txt','w')
        handle_out.write("ABI trace data\nAndrew's text out version (sqrt(5) + 1)/2\n\n")

        # print metadata header
        for key in ['data collection finish date','run start time', 'data collection start date', 'polymer', 'data collection start time', 'well', 'run finish time', 'data collection finish time', 'baseorder', 'run start date', 'dye', 'model', 'run finish date']:
            handle_out.write(key + ': ' + str(test.data[key]) +'\n')

        trace_data1 = test.tags['DATA9'].tag_data
        trace_data2 = test.tags['DATA10'].tag_data
        trace_data3 = test.tags['DATA11'].tag_data
        trace_data4 = test.tags['DATA12'].tag_data

        baseorder = test.data['baseorder']
        trace_len = len(test.tags['DATA9'].tag_data)

        j=0
        seq_len = len(test.seq)
        next_base=test.seq[j]

        handle_out.write('\n')
        handle_out.write('\tSeq\t\tTraces\t\t\t\t\tSequence\n')
        handle_out.write('\t\t')
        for base in baseorder:
            handle_out.write(base + '\t')
        handle_out.write('\n')

        # write out tab delimeted sequence
        for i in range(trace_len):
            handle_out.write(str(i) + '\t')
            if i in test.data['tracepeaks']:
                handle_out.write(next_base + '\t')
                j += 1
                if j<seq_len:
                    next_base = test.seq[j]
                else:
                    next_base = 'N'
            else:
                handle_out.write('.\t')
            handle_out.write(str(trace_data1[i]) + '\t')
            handle_out.write(str(trace_data2[i]) + '\t')
            handle_out.write(str(trace_data3[i]) + '\t')
            if i>=len(test.seq):
                handle_out.write(str(trace_data4[i]) + '\n')
            else:
                handle_out.write(str(trace_data4[i]) + '\t\t')
                handle_out.write(str(i+1) + '\t')
                handle_out.write(test.seq[i])
                handle_out.write('\n')

        # write fasta files into same directory for handling by Igor's script
        print  directory + root_name + '.txt'
        fasta_out = open( directory + root_name + '.fasta','w')
        abi_in = open(filename,'rb')
        sequences = SeqIO.parse(abi_in,'abi')
        SeqIO.write(sequences, fasta_out,'fasta')
        fasta_out.close()
        abi_in.close()
