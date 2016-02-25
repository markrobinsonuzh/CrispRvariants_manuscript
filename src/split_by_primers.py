import regex
import os, sys
import gzip
from argparse import ArgumentParser


def sq_to_rc(seq):
  seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
  return "".join([seq_dict[base] for base in reversed(seq)])


def next_record(fq):
  first_id = fq.readline()
  if (first_id == ""):
    return(False)
  seq = fq.readline()
  second_id = fq.readline()
  quals = fq.readline()
  return((first_id, seq, second_id, quals))
  
  
def match_primer(record, res, rc = False):
  sq = record[1].strip()
  if (rc == True): 
    sq = sq_to_rc(record[1])
  for i, re_exp in enumerate(res):
    mm = regex.match(re_exp, sq)
    if mm:
      return(i)
  return(-1)
   
   
def setup_outdirs(out_dir, primer_names):
  out_dirs = [os.path.join(out_dir, nm) for nm in primer_names]
  for out_dir in out_dirs:
    if not os.path.exists(out_dir):
      os.mkdir(out_dir)
  return(out_dirs)
  
  
def print_tallies(primer_names, tallies):
  for (nm, t) in zip(primer_names, tallies):
    print(nm, t)
  print("not found: %s" % tallies[-1])
  return()
  
     
def close_all(files):
  for f in files:
    f.close()

   
def main_paired(f1_fname, f2_fname, out_dir, fwds, primer_names):
  f1 = gzip.open(f1_fname, "r")
  f2 = gzip.open(f2_fname, "r")
  
  r1 = next_record(f1)
  r2 = next_record(f2)
  
  fwds = [fwd for (id, fwd) in primers]
  primer_names = [id for (id, fwd) in primers]
  
  fwd_only = [regex.compile("(%s){e<=1}" % fwd) for fwd in fwds]

  out_dirs = setup_outdirs(out_dir, primer_names)

  out_fwd = [gzip.open(os.path.join(od, os.path.basename(f1_fname)), "wb") \
              for od in out_dirs]
  out_rev = [gzip.open(os.path.join(od, os.path.basename(f2_fname)), "wb") \
              for od in out_dirs]

  i = 0 
  done = False
  tallies = [0] * (len(primer_names) + 1)  

  while (done == False):
    re1 = match_primer(r1, fwd_only)
    if (re1 != -1):
      tallies[re1] += 1
      out_fwd[re1].writelines(r1)
      out_rev[re1].writelines(r2)
    else:
      tallies[-1] += 1    
        
    r1 = next_record(f1)
    r2 = next_record(f2)
    
    i += 1
    if i % 10000 == 0:
      print(i, tallies)
      
    if (r1 == False):
      done = True
        
  fs = out_fwd + out_rev + [f1, f2]
  close_all(fs)
  
  print_tallies(primer_names, tallies)  
  

def main_single(f1_fname, out_dir, fwds, primer_names):
  f1 = gzip.open(f1_fname, "r")
  r1 = next_record(f1)
  
  fwd_only = [regex.compile("(%s){e<=1}" % fwd) for fwd in fwds]

  out_dirs = setup_outdirs(out_dir, primer_names)
  out = [gzip.open(os.path.join(od, os.path.basename(f1_fname)), "wb") \
              for od in out_dirs]

  i = 0 
  done = False
  tallies = [0] * (len(primer_names) + 1)  

  while (done == False):
    re1 = match_primer(r1, fwd_only)
    if (re1 != -1):
      tallies[re1] += 1
      out[re1].writelines(r1)
    else:
      tallies[-1] += 1    
        
    r1 = next_record(f1)
    
    i += 1
    if i % 10000 == 0:
      print(i, tallies)
      
    if (r1 == False):
      done = True
        
  close_all(out)
  f1.close()
  
  print_tallies(primer_names, tallies)

def main_compare():
  f = open("../annotation/Shah_metadata_edited.txt", "r") 
  header = f.next()
  
  # Strip adapters, reverse complement the reverse primers  
  primers = [(l[0], l[3], sq_to_rc(l[4])) \
     for l in (line.strip().split("\t") for line in f)]
   
  # Test different primer matching rules
  strict = [regex.compile("(%s){e<=1}.*(%s){e<=1}" % (fwd, rev)) for \
         (id, fwd, rev) in primers]
  
  tolerant = [regex.compile(".*(%s){e<=2}.*(%s){e<=2}.*" % (fwd, rev)) for \
         (id, fwd, rev) in primers]
  
  
  # Setup directories and outfiles
  dirs = [os.path.join("../merged_split", x) for x in \
          ["strict_pear", "tolerant_pear"]]
  
  for dr in dirs:
    if not os.path.exists(dr):
      os.mkdir(dr)
  
  out_f_sets = [[open(os.path.join(dr, "%s.fastq" % id), "w") \
                 for (id, fwd, rev) in primers] for dr in dirs]

  fq = gzip.open("../fastq/SRR1769728_merged_pear.assembled.fastq.gz", "r")
  
  
  # Match primers
  r1 = next_record(fq)

  done = False
  while (done == False):
    for (res, out) in zip((strict, tolerant), out_f_sets):
      re1 = match_primer(r1, res)
      if (re1 != -1):
        out[re1].writelines(r1)
          
    r1 = next_record(fq)
    if (r1 == False):
      done = True
      close_all(out)
  fq.close()

  
    
if __name__ == "__main__":  

  parser = ArgumentParser()
  parser.add_argument('-r1', dest='r1', default="", help='forward reads')
  parser.add_argument('-r2', dest='r2', default="", help='reverse reads')
  parser.add_argument('-o', dest='out_dir', default="", help='output directory')
  parser.add_argument('-p', dest='primers', default="", 
                      help='primers filename (id fwd, tab separated)')
  parser.add_argument('-c', dest='compare', action = "store_true", default=False, 
                      help='run the primer comparison used in CrispRVariants paper')
  
  args = parser.parse_args()
  if (args.compare):
    print("Running primer comparison")
    main_compare()
    sys.exit()

  primers = [l.strip().split("\t") for l in open(args.primers, "r")]
  fwds = [fwd for (id, fwd) in primers]
  primer_names = [id for (id, fwd) in primers]
  
  if (args.r2 != ""):
    main_paired(args.r1, args.r2, args.out_dir, fwds, primer_names)
  else:
    main_single(args.r1, args.out_dir, fwds, primer_names)  
