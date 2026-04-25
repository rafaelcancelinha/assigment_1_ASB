#!/usr/bin/env python3
"""
Complete phylogenetic pipeline:
  [Full mode]     python pipeline.py <database> <accessions.txt>
       Fetch → Clean → Align → ModelTest → RAxML → MrBayes

  [Aligned mode]  python pipeline.py <concatenated.fasta>
       ModelTest → RAxML → MrBayes  (starts from aligned FASTA)
"""
import subprocess
import os
import sys
import time
import re
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment


def extract_model(filename):
    """Reads the best model from ModelTest-NG."""
    if not os.path.exists(filename):
        print(f"Warning: {filename} not found. Using GTR+G.")
        return "GTR+G"
    try:
        with open(filename) as f:
            for line in f:
                if "Best-fit model:" in line or "Model:" in line:
                    model = line.split(":")[-1].strip()
                    # RAxML-NG prefers GTR if the model is TIM or TVM
                    if model.startswith(("TIM", "TVM")):
                        return "GTR+G"
                    return model
    except Exception:
        pass
    return "GTR+G"


def sanitize_nexus_name(name):
    """Removes problematic characters for MrBayes."""
    name = name.replace("'", "").replace('"', "").replace(" ", "_")
    name = re.sub(r'[^A-Za-z0-9_.\-]', '_', name)
    return name


def convert_to_nexus(fasta_input, nexus_output):
    """Converts aligned FASTA to NEXUS format compatible with MrBayes."""
    records = []
    for rec in SeqIO.parse(fasta_input, "fasta"):
        rec.id          = sanitize_nexus_name(rec.id)
        rec.name        = rec.id
        rec.description = ""
        rec.annotations["molecule_type"] = "DNA"
        records.append(rec)
    AlignIO.write(MultipleSeqAlignment(records), nexus_output, "nexus")

    with open(nexus_output, 'r') as f:
        content = f.read()
    content = content.replace("'", "").replace('"', "")
    with open(nexus_output, 'w') as f:
        f.write(content)


def run_mrbayes(nexus_file, model):
    """Runs MrBayes with automatic stopping rule."""
    print("\n--- Step: MrBayes (Optimized with Auto-stop) ---")
    model_mb = model.split("+")[0]
    nst = "6" if model_mb in ["GTR", "TIM", "TIM1", "TIM2", "TIM3", "TVM", "TrN", "TrNef", "TVMef", "SYM"] else "2"
    mrbayes_file = nexus_file.replace(".nex", ".mrbayes.nex")

    block = f"""begin mrbayes;
  set autoclose=yes;
  lset nst={nst} rates=gamma;
  mcmc ngen=10000000 nchains=4 nruns=2 samplefreq=1000 burninfrac=0.25 stoprule=yes stopval=0.01;
  sumt burninfrac=0.25 contype=halfcompat;
  quit;
end;"""

    with open(nexus_file, 'r') as f:
        content = f.read().rsplit("end;", 1)[0]
    with open(mrbayes_file, "w") as f:
        f.write(content + "\nend;\n" + block)

    subprocess.run(["mb", mrbayes_file], check=True)
    return mrbayes_file + ".con.tre"


def run_raxml(aligned_fasta, model, prefix):
    """Runs RAxML-NG with automatic bootstraps (autoMRE) and 8 threads."""
    print("\n--- Step: RAxML-NG (autoMRE, 8 Threads) ---")
    subprocess.run([
        "raxml-ng",
        "--all",
        "--msa",      aligned_fasta,
        "--model",    model,
        "--bs-trees", "autoMRE",
        "--threads",  "8",
        "--seed",     "12345",
        "--prefix",   prefix,
        "--redo"
    ], check=True)
    return f"{prefix}.raxml.support"


def main():
    n_args = len(sys.argv) - 1

    # ------------------------------------------------------------------ #
    #  ALIGNED MODE – 1 argument: <concatenated.fasta>                   #
    # ------------------------------------------------------------------ #
    if n_args == 1:
        ALIGNED_INPUT = sys.argv[1]
        base    = os.path.splitext(os.path.basename(ALIGNED_INPUT))[0]
        workdir = f"tree_{base}"
        os.makedirs(workdir, exist_ok=True)
        NEXUS = os.path.join(workdir, f"{base}.nex")

        # Convert to NEXUS
        convert_to_nexus(ALIGNED_INPUT, NEXUS)

        # ModelTest-NG
        print("\n--- Step: ModelTest-NG (12 Threads) ---")
        subprocess.run(
            ["modeltest-ng", "-i", ALIGNED_INPUT, "-d", "nt", "-p", "12"],
            check=True
        )
        model = extract_model(f"{ALIGNED_INPUT}.out")
        print(f"\n  Selected model: {model}")

        # RAxML-NG
        ml_tree_path = run_raxml(ALIGNED_INPUT, model, os.path.join(workdir, "raxml"))

        # MrBayes
        bayes_tree_path = run_mrbayes(NEXUS, model)

        print("\n" + "="*60)
        print("    Pipeline completed!")
        print(f"   Results folder      : {workdir}")
        print(f"   ML Tree (RAxML)     : {ml_tree_path}")
        print(f"   Bayesian Tree       : {bayes_tree_path}")
        print("="*60)

    # ------------------------------------------------------------------ #
    #  FULL MODE – 2 arguments: <database> <accessions.txt>              #
    # ------------------------------------------------------------------ #
    elif n_args == 2:
        DATABASE, ACCESSIONS_FILE = sys.argv[1], sys.argv[2]

        base    = os.path.splitext(os.path.basename(ACCESSIONS_FILE))[0]
        workdir = f"tree_{base}"
        os.makedirs(workdir, exist_ok=True)

        RAW     = os.path.join(workdir, "temp.fasta")
        CLEAN   = os.path.join(workdir, "clean.fasta")
        ALIGNED = os.path.join(workdir, "aligned.fasta")
        NEXUS   = os.path.join(workdir, "aligned.nex")

        # 1. Fetch sequences
        from fetch_sequences import fetch_all
        fetch_all(DATABASE, ACCESSIONS_FILE, RAW)

        # 2. Cleaning
        from simple_fasta import simplify_fasta
        simplify_fasta(RAW, CLEAN)

        # 3. Alignment
        with open(ALIGNED, "w") as out:
            subprocess.run(["mafft", "--auto", CLEAN], stdout=out, check=True)

        # 4. Convert to NEXUS
        convert_to_nexus(ALIGNED, NEXUS)

        # 5. ModelTest-NG
        print("\n--- Step: ModelTest-NG (12 Threads) ---")
        subprocess.run(
            ["modeltest-ng", "-i", ALIGNED, "-d", "nt", "-p", "12"],
            check=True
        )
        model = extract_model(f"{ALIGNED}.out")
        print(f"\n  Selected model: {model}")

        # 6. RAxML-NG
        ml_tree_path = run_raxml(ALIGNED, model, os.path.join(workdir, "raxml"))

        # 7. MrBayes
        bayes_tree_path = run_mrbayes(NEXUS, model)

        print("\n" + "="*60)
        print("    Pipeline completed!")
        print(f"   Results folder      : {workdir}")
        print(f"   ML Tree (RAxML)     : {ml_tree_path}")
        print(f"   Bayesian Tree       : {bayes_tree_path}")
        print("="*60)

    # ------------------------------------------------------------------ #
    #  Invalid arguments                                                 #
    # ------------------------------------------------------------------ #
    else:
        print("Usage:")
        print("  python pipeline.py <concatenated.fasta>          # Already aligned FASTA")
        print("  python pipeline.py <database> <accessions.txt>   # Full pipeline")
        sys.exit(1)


if __name__ == "__main__":
    main()

