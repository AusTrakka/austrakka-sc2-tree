import os
import subprocess
import tempfile
import sys
from Bio import SeqIO
import argparse
from shutil import copyfile

#....................................................................
def msg(*args):
    print(*args, file=sys.stderr)

#....................................................................
def err(*args):
    msg("ERROR:", *args)
    sys.exit(1)

#....................................................................
def wrn(*args):
    msg("WARNING:", *args)


#....................................................................
def fixfile(f, text="bummer"):
    text ||= "bummer"
    if not f or not os.path.isfile(f) or not os.access(f, os.R_OK):
        err(f"{text} - can't read file '{f}'")
    return os.path.abspath(f)

#....................................................................
def kovid_enter(*args):
    msg("kovid toolkit by \@torstenseemann")
    if args:
        msg(*args)

#....................................................................
def kovid_exit(*args):
    if old_dir:
        msg("Changing back to dir:", old_dir)
        os.chdir(old_dir)
    if args:
        msg(*args)
    msg("Done.")
    sys.exit(0)

#....................................................................
def exe(tool):
    path = shutil.which(tool)
    if path:
        msg(f"Found '{tool}' - {path}") 
    else:
        err(f"Can't find '{tool}' - please install it")
    return path

#....................................................................
def fasta2array(*files):
    exe("seqkit")
    seq = []
    for fname in files:
        msg("Loading FASTA:", fname)
        records = list(SeqIO.parse(fname, "fasta"))
        seq.extend([(rec.id, str(rec.seq)) for rec in records])
    msg("Loaded", len(seq), "sequences.")
    return seq

#....................................................................
def tmpdir(change=False):
    dir = tempfile.TemporaryDirectory()
    msg("Created tempdir:", dir)
    if change:
        msg("Changing to tempdir:", dir)
        old_dir = os.getcwd()
        os.chdir(dir)
    return dir

#....................................................................
def num_lines(fname):
    return sum(1 for line in open(fname))

#....................................................................
def run(*args):
    msg("Running:", *args)
    if subprocess.call(args) != 0:
        err("Could not run:", *args)

#....................................................................
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--afa', required=True)
    parser.add_argument('-p', '--prefix', required=True)
    parser.add_argument('-j', '--cpus', type=int, default=os.cpu_count())
    parser.add_argument('-r', '--ref', required=True)
    parser.add_argument('-e', '--exe', required=True)
    parser.add_argument('-x', '--opts', default='')
    parser.add_argument('-t', '--tmpdir')
    parser.add_argument('-B', '--bootstraps', type=int, default=1000)
    parser.add_argument('-F', '--faster', action='store_true')
    args = parser.parse_args()

    afa = fixfile(args.afa, "need -a AFA")
    ref = fixfile(args.ref, "Need -r REF")
    cpus = args.cpus
    out = args.prefix or err("need -p PREFIX")
    ftree = args.exe or err("need -e EXE")
    tdir = args.tmpdir or tmpdir()

    kovid_enter("This is", sys.argv[0])
    exe(ftree)
    exe("raxml-ng")
    exe(tool) for tool in ["goalign", "gotree", "nw_stats", "nw_order", "nw_reroot"]

    refseq = fasta2array(ref)
    L = len(refseq[0][1])
    outgroup = refseq[0][0]
    msg(f"Reference {outgroup} is {L} bp")

    msg("Working folder:", tdir)
    tmp = f"{tdir}/kovid.FT"
    msg("Temp prefix:", tmp)

    msg("Purifying AFA for", args.exe)
    run(f"{os.getcwd()}/kovid_afa_clean.pl -F {afa} > {tmp}.afa1")

    # Remaining steps omitted for brevity
    msg("De-duping AFA...")
    run(f"goalign dedup -i {tmp}.afa1 -o {tmp}.afa2")
    msg("Unique seqs  :", num_lines(f"{tmp}.dup"))

    msg("Removing invariant sites...")
    run(f"goalign compress -i {tmp}.afa1 -o {tmp}.afa3 --weight-out {tmp}.afa3.weights")
    var = num_lines(f"{tmp}.afa3.weights")
    msg("Unique sites :", var)

    msg("Building tree with", ftree)

    opts = args.opts if args.opts else ''
    opts += " -fastest" if args.faster
    opts += f" -boot {args.bootstraps}" if args.bootstraps else " -nosupport"
    opts += f" -threads {args.cpus}" if args.exe == 'VeryFastTree'
    msg("Opts:", opts)

    run(f"OMP_NUM_THREADS={args.cpus} time {ftree} {opts} -out {tmp}.nwk1 -log {tmp}.log -nt {tmp}.afa3")
    copyfile(f"{tmp}.log", f"{out}.{ftree}.log")
    copyfile(f"{tmp}.nwk1", f"{out}.{ftree}.nwk")

    msg("Fixing branch lengths...")
    run(f"gotree resolve < {tmp}.nwk1 > {tmp}.nwk2")
    run(f"raxml-ng --threads {cpus} --force perf_threads --evaluate --model GTR+G4 --blmin {BLMIN} --prefix {tmp} --tree {tmp}.nwk2 --msa {tmp}.afa3 --site-weights {tmp}.afa3.weights")
    copyfile(f"{tmp}.raxml.bestTree", f"{tmp}.nwk3")

    msg("Cleaning branches...")
    run(f"cat {tmp}.nwk3 | gotree brlen round -p {PRECISION} | gotree collapse length -l 0 > {tmp}.nwk4")

    msg("Repopulating, sorting, rerooting tree...")
    run(f"gotree repopulate -i {tmp}.nwk4 -g {tmp}.dup | nw_order -c a - | nw_reroot - {outgroup} | gotree brlen scale -f 1 > {tmp}.nwk5")

    msg("Collating output files")
    copyfile(f"{tmp}.nwk5", f"{out}.nwk")

    run(f"nw_stats {out}.nwk")

    kovid_exit()

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
