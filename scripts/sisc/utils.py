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
