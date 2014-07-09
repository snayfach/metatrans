# Libraries
import os, optparse, gzip, subprocess
import datetime as dt

# Parse command line args
parser          = optparse.OptionParser(usage = "Usage: python metatrans.py [-options] <p_reads> <p_orfs>")
parser.add_option("-m", dest="method", default='Prodigal', help="translation method [6FT, 6FT-split, Prodigal] Default = Prodigal")
parser.add_option("-l", dest="orf_len", default='1', help="minimum orf length to keep. Default = 1")
parser.add_option("-t", dest="timeit", default=False, action='store_true', help="time execution. Default = False")
try:
    (options, args) = parser.parse_args()
    p_reads         = args[0]
    p_orfs          = args[1]
    method          = options.method
    orf_len         = int(options.orf_len)
    timeit          = options.timeit
#    p_reads         = "./example/test.fna.gz" # example params
#    p_orfs          = "./example/test.faa.gz" # example params
#    method          = "6FT" # example params
#    orf_len         = 15 # example params
except Exception:
    print "\nIncorrect options/parameters!"
    print "Usage: python metatrans.py [-options] <p_reads> <p_orfs>"
    print ""
    quit()
if method not in ['6FT', '6FT-split', 'Prodigal']:
    print 'Error: Incorrect value for -m option. Must be one of: 6FT, 6FT-split, Prodigal'
    exit()
if orf_len > 1 and method != '6FT-split':
    print 'Warning: ORF length filtering only valid when using -m 6FT-split.'

# Functions
def transeq(p_reads):
    """
        Run transeq on p_reads; Return results in pipe
    """
    gz_in = True if p_reads.split('.')[-1] == 'gz' else False
    if gz_in:
        command = "zcat %s | ./transeq -trim -frame=6 -sformat1 pearson -osformat2 pearson stdin stdout" % p_reads
    else:
        command = "zcat  | ./transeq -trim -frame=6 -sformat1 pearson -osformat2 pearson %s stdout" % p_reads
    pipe = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=open('/dev/null', 'w'))
    return pipe.stdout

def six_frame_trans(p_reads, p_orfs, method, orf_len):
    pipe = transeq(p_reads)
    f_orfs = gzip.open(p_orfs, 'w') if p_orfs.split('.')[-1] == 'gz' else open(p_orfs, 'w')
    if method == '6FT':
        write_no_split(pipe, f_orfs)
    elif method == '6FT-split':
        write_no_split(pipe, f_orfs, orf_len)

def write_no_split(pipe, f_out):
    """
        Write sequences from pipe to f_out
    """
    id = ''
    seq = ''
    frame = 0
    # init seq id
    for line in pipe:
        id = line.rstrip().split()[0]
        break
    # loop over lines
    for line in pipe:
        if line[0] == '>':
            f_out.write(id+'_'+str(frame)+'\n'+seq+'\n')
            id = line.rstrip().split()[0]
            seq = ''
        else:
            seq += line.rstrip()
    f_out.write(id+'_'+str(frame)+'\n'+seq+'\n')

def write_with_split(pipe, f_out, orf_len):
    """
        Split sequences from pipe on stop codons; write split seqs to f_out if they meet or exceed orf_len
    """
    id = ''
    seq = ''
    # init seq id
    for line in pipe:
        id = line.rstrip().split()[0]
        break
    # loop over lines
    for line in pipe.stdout:
        if line[0] == '>':
            for frame, split_seq in enumerate(seq.split('*')):
                if len(split_seq) >= orf_len:
                    f_out.write(id+'_'+str(frame)+'\n'+split_seq+'\n')
            id = line.rstrip().split()[0]
            seq = ''
        else:
            seq += line.rstrip()
    # write last seq
    f_out.write(id+'_'+str(frame)+'\n'+seq+'\n')


def run_prodigal(p_reads, p_orfs):
    """
        Run prodigal on p_reads. Write results to p_orfs
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
start = dt.datetime.now()
if method in ['6FT','6FT-split']:
    six_frame_trans(p_reads, p_orfs, method, orf_len)
elif method == 'Prodigal':
    run_prodigal(p_reads, p_orfs)
stop = dt.datetime.now()
if timeit:
    print 'Time elapsed (s):', (stop - start).seconds




