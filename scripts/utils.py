def run_test(kernel, ts, nx, ny, nz, nt, exe_cmd, outfile, target_dir, is_dp=1, dry_run=0, **kwargs):
  import os
  import subprocess
  from string import Template
  from scripts.utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests $ntests --disable-source-point --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz  --verbose $verbose --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs --mwd-type $mwdt --bsx $bsx --num-wavefronts $nwf --verify $verify | tee $outfile""")

  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

  defaults = {'is_dp':1, 'tgs':1, 'cs':8192, 'mwdt':1, 'npx':1, 'npy':1, 'npz':1, 'nwf':1,
              'bsx':1e6, 'ntests':2, 'alignment':16, 'verify':0, 'verbose':1}
  defaults.update(kwargs)
  job_cmd = job_template.substitute(nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, 
                                    exec_path=exec_path, target_dir=target_dir, **defaults)


  job_cmd = exe_cmd + job_cmd
 
  print job_cmd
  if(dry_run==0): sts = subprocess.call(job_cmd, shell=True)

  return job_cmd

def select_fields(data, rows=[], cols=[]):

  sel_data = []
  # select rows
  if rows:
    for k in data:
      for row in rows:
        cmp_res = [k[field] == val for field, val in row.iteritems()]
        if all(cmp_res):
          if cols:
            tup = {}
            for col in cols:
              tup[col] = k[col]
            sel_data.append(tup)
          else: # select all columns
            sel_data.append(dict(k))
        

  return sel_data


def create_project_tarball(dest_dir, fname):
  import tarfile, glob
  import os

  nl = ["src/kernels/*.c", "src/*.c", "src/*.h", "scripts/*.py", "scripts/*/*.py", "make.inc", "Makefile"]
  nl = [glob.glob(n) for n in nl]
  nl = [n for nn in nl for n in nn]

  out_dir = os.path.join(os.path.abspath("."),dest_dir)
  ensure_dir(out_dir)
  out_name = os.path.join(out_dir, fname+".tar.gz")

  print "Writing project files to:" + out_name
  with tarfile.open(out_name, "w:gz") as tar:
    for n in nl:
      print "Adding to the tar file: " + n
      tar.add(n)


def ensure_dir(d):
  import os, errno
  try:
    os.makedirs(d)
  except OSError as exc:
    if exc.errno == errno.EEXIST:
      pass
    else: raise

def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data


