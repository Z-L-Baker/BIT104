#!/bin/bash
#SBATCH --partition=defq       # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=2      #
#SBATCH --mem=2000     # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err         # redirect stderr to this file
#SBATCH --output=%J.out        # redirect stdout to this file

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

# Write jobscript to output file (good for reproducability)
cat $0

workdir=$(pwd)

awk -F'\t' -v OFS='\t' \
    -v HUMAN_TAG="Hsap" \
    -v STRONG_HUMAN=1e-20 \
    -v MAX_HUMAN_E=1e-10 \
    -v MAX_NONHUMAN_E=1e-10 '
{
  gene=$1; hit=$2; e=$3+0; sp=$5;

  # remove the gene .t suffix
  sub(/\.t[0-9]+$/, "", gene)

  # --- Extract protein symbol ---
  prot=hit
  sub(/.*\|/,"",prot)       # remove everything before last "|"
  sub(/_[A-Z]+$/,"",prot)   # remove species suffix (e.g. _HUMAN)

  # --- Extract UniProt ID ---
  split(hit, parts, "|")
  uniprot = parts[2]

  genes[gene]=1;

  # Store best human vs non-human (protein symbol)
  if (sp==HUMAN_TAG) {
    if (!(gene in he) || e < he[gene]) { he[gene]=e; hh[gene]=prot; hhu[gene]=uniprot }
  } else {
    if (!(gene in nhe) || e < nhe[gene]) { nhe[gene]=e; nhh[gene]=prot; nhhu[gene]=uniprot }
  }
}
END {
  # Print headers
  print "gene_id","hit","evalue","decision" > "gene2protID.tsv"
  print "gene_id","uniprot","evalue","decision" > "gene2uniprotID.tsv"

  for (g in genes) {
    Hprot=(g in hh)?hh[g]:"";
    He=(g in he)?he[g]:"";
    Huni=(g in hhu)?hhu[g]:"";

    NHprot=(g in nhh)?nhh[g]:"";
    NHe=(g in nhe)?nhe[g]:"";
    NHuni=(g in nhhu)?nhhu[g]:"";

    decision="";
    bestProt=""; bestUni=""; bestE="";

    if (He != "") {
      # --- case: we have a human hit ---
      if (Hprot == NHprot) {
        decision="same protein → human kept"
        bestProt=Hprot; bestUni=Huni; bestE=He
      }
      else if (He <= STRONG_HUMAN) {
        decision="strong human"
        bestProt=Hprot; bestUni=Huni; bestE=He
      }
      else if (He <= MAX_HUMAN_E) {
        decision="human ok"
        bestProt=Hprot; bestUni=Huni; bestE=He
      }
      else {
  	# human weaker than ok, compare against non-human
  	if (NHprot != "" && NHe <= MAX_NONHUMAN_E && NHe < He) {
	  decision="human weak, better non-human used"
    	  bestProt=NHprot; bestUni=NHuni; bestE=NHe
  	} else {
    	  decision="human weak, no acceptable hit"
          bestProt=""; bestUni=""; bestE=""
  	}
      }
    }
    else if (NHprot != "" && NHe <= MAX_NONHUMAN_E) {
      # --- no human hit, but good non-human ---
      decision="no human, non-human used"
      bestProt=NHprot; bestUni=NHuni; bestE=NHe
    }
    else {
      # --- nothing useful ---
      decision="no hits"
      bestProt=""; bestUni=""; bestE=""
    }

    # Print to protein symbol file
    printf "%s\t%s\t%s\t%s\n", g, bestProt, bestE, decision >> "gene2protID.tsv"

    # Print to UniProt ID file
    printf "%s\t%s\t%s\t%s\n", g, bestUni, bestE, decision >> "gene2uniprotID.tsv"
  }
}
' blast_best.tsv
sort -k1,1 gene2protID.tsv -o gene2protID.tsv
sort -k1,1 gene2uniprotID.tsv -o gene2uniprotID.tsv
