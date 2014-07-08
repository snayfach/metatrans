# Libraries
import os, optparse, gzip, subprocess
import Bio.SeqIO
import Bio.Seq
import Bio.SeqUtils

# Parse command line args
parser          = optparse.OptionParser(usage = "Usage: python metatrans.py [-options] <p_reads> <p_orfs>")
parser.add_option("-m", dest="method", default='6FT', help="translation method [6FT, 6FT-split, Prodigal] Default = 6FT")
parser.add_option("-l", dest="orf_len", default='1', help="minimum orf length to keep. Default = 1")

try:
    (options, args) = parser.parse_args()
    p_reads         = args[0]
    p_orfs          = args[1]
    method          = options.method
    orf_len         = int(options.orf_len)
#    p_reads         = "./example/test.fna.gz" # example params
#    p_orfs          = "./example/test.faa.gz" # example params
#    method          = "6FT-split" # example params
#    orf_len         = 1 # example params
except Exception:
    print "\nIncorrect options/parameters!"
    print "Usage: python metatrans.py [-options] <p_reads> <p_orfs>"
    print ""
    quit()

# Functions
def six_frame_trans(p_reads):
    """
        Translate read into all 6 possible reading frames
    """
    f_in = gzip.open(p_reads) if p_reads.split('.')[-1] == 'gz' else open(p_reads)
    for record in Bio.SeqIO.parse(f_in, 'fasta'):
        read_id  = record.id
        orf_id   = 0
        frame_id = 0
        my_seq   = record.seq
        rev_comp = Bio.Seq.reverse_complement(my_seq)
        seq_len  = len(record.seq)
        for start in range(3):
            orf_id += 1
            seq_id  = read_id+'_'+str(orf_id)+'_'+str(frame_id)
            stop    = seq_len - ((seq_len - start) % 3)
            orf     = Bio.Seq.translate(my_seq[start:stop])
            yield (seq_id, str(orf))
        for start in range(3):
            orf_id += 1
            seq_id  = read_id+'_'+str(orf_id)+'_'+str(frame_id)
            stop = seq_len - ((seq_len - start) % 3)
            Bio.Seq.translate(rev_comp[start:stop])
            yield (seq_id, str(orf))

def split_orf(seq):
    """
        Translate read into all 6 possible reading frames
    """
    orf_id   = seq[0]
    orf_seq  = seq[1]
    frame_id = 0
    for orf_split in orf_seq.split('*'):
        frame_id += 1
        seq_id  = orf_id.rstrip('0')+str(frame_id)
        yield (seq_id, orf_split)

def run_prodigal(p_reads, p_orfs):
    """
        Run prodigal
    """
    p_prodigal = './Prodigal-2.60/prodigal' # path to binary
    p_tmp = p_orfs + 'tmp' # tmp file
    # run prodigal
    gz_in = True if p_reads.split('.')[-1] == 'gz' else False
    if gz_in: command = 'zcat %s | %s -a %s -p meta -q -o /dev/null' % (p_reads, p_prodigal, p_tmp)
    else: command = '%s -i %s -a %s -p meta -q -o /dev/null' % (p_prodigal, p_reads, p_tmp)
    process = subprocess.Popen(command, shell=True)
    process.wait()
    # process output file
    gz_out = True if p_orfs.split('.')[-1] == 'gz' else False
    if gz_out:
        f_in = open(p_tmp)
        f_out = gzip.open(p_orfs, 'w')
        for line in f_in:
            f_out.write(line)
        f_in.close()
        f_out.close()
        os.remove(p_tmp)
    else:
        os.rename(p_tmp, p_orfs)

# Main
if method == '6FT':
    f_out = gzip.open(p_orfs, 'w') if p_orfs.split('.')[-1] == 'gz' else open(p_orfs, 'w')
    for orf in six_frame_trans(p_reads):
        f_out.write('>'+orf[0]+'\n')
        f_out.write(orf[1]+'\n')
elif method == '6FT-split':
    f_out = gzip.open(p_orfs, 'w') if p_orfs.split('.')[-1] == 'gz' else open(p_orfs, 'w')
    for orf in six_frame_trans(p_reads):
        for frame in split_orf(orf):
            if len(frame[1]) >= orf_len:
                f_out.write('>'+frame[0]+'\n')
                f_out.write(frame[1]+'\n')
elif method == 'Prodigal':
    run_prodigal(p_reads, p_orfs)




