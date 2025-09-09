#!/bin/bash
#SBATCH --job-name=merge
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
blastdir=${workdir}/rawreads/transpiper/outdir/blast_results
upidir=${workdir}/rawreads/transpiper/outdir/annotations/upimapi
iprdir=${workdir}/interproscan
endir=${workdir}/eggnog

## Make blast file

echo -e "gene\tsubject\tevalue\tscore\tspecies" > blast_all.tmp

for f in ${blastdir}/e_crypt_gen_guided_*.tsv
	do

	R1=$(basename "$f" | cut -f1 -d.)
	base=$(echo "$R1" | sed 's/e_crypt_gen_guided_//')

	species=$(echo "$base" | sed 's/_blp//')

	awk -v sp="$species" 'BEGIN{OFS="\t"} {print $1,$2,$11,$12,sp}' "$f"
done >> blast_all.tmp

awk 'NR==1 {print; next} {
    key=$1 FS $5
    if (!(key in best) || $3 < eval[key] || ($3 == eval[key] && $4 > score[key])) {
        best[key]=$0
        eval[key]=$3
        score[key]=$4
    }
}
END {
    print "gene\tsubject\tevalue\tscore\tspecies"
    for (k in best) print best[k]
}' blast_all.tmp | sort -k1,1 > blast_best.tsv

## Get GO terms for blast results from upimapi results

for upi in ${upidir}/*/e_crypt_gen_guided_*_upimapi.tsv; do
    awk -F'\t' '
    NR==1 {
        # Find the column numbers by header names
        for (i=1;i<=NF;i++) {
            if ($i=="Entry") entry_col=i
            if ($i=="Gene Ontology (GO)") go_col=i
        }
        next
    }
    $entry_col != "" && $go_col != "" { print $entry_col"\t"$go_col }
    ' "$upi"
done | sort -u > all_entry_go.tsv

awk -F'\t' '
NR==FNR {
    go[$1]=$2
    next
}
NR==1 {
    print $0"\tGO_terms"; next
}
{
    split($2,a,"|")   # split subject by "|"
    uniprot=a[2]      # <-- this is the correct UniProt ID
    print $0"\t"(uniprot in go ? go[uniprot] : "")
}' OFS='\t' all_entry_go.tsv blast_best.tsv > blast_with_go.tsv


# Make a simplified file with only ID and GO
awk '
{
    id = $1;
    if (!(id in seen)) {
        order[++n] = id;   # keep track of first appearance
        seen[id] = 1;
    }
    # extract GO terms
    while (match($0, /\[GO:[0-9]+\]/)) {
        term = substr($0, RSTART+1, RLENGTH-2);
        key = id SUBSEP term;   # unique key per ID + GO term
        if (!(key in goterm)) {
            goterm[key] = 1;
            terms[id] = terms[id] term ",";
        }
        $0 = substr($0, RSTART+RLENGTH);
    }
}
END {
    for (i=1; i<=n; i++) {
        id = order[i];
        out = terms[id];
        sub(/,$/, "", out);
        print id "\t" out;
    }
}' blast_with_go.tsv > blast_go_collapsed.tsv


## Make interproscan file

awk -F"\t" '
{
    id = $1
    if ($14 != "-") {
        n = split($14, arr, "|")
        for (i=1; i<=n; i++) {
            # remove "(source)" e.g. GO:0005515(InterPro)
            goid = arr[i]
            sub(/\(.*\)$/, "", goid)
            if (!seen[id, goid]++) {
                goterms[id] = (goterms[id] ? goterms[id] "," goid : goid)
            }
        }
    }
}
END {
    OFS="\t"
    for (id in goterms) {
        print id, goterms[id]
    }
}' ${iprdir}/E_crypticus_protein.fasta_1.tsv | sort -k1,1 > ipr_collapsed.tsv




## Make eggnog file

grep -v '^##' ${endir}/ec_eggnog.emapper.annotations > ec_eggnog.tsv

# Extract only the header and data lines, then print query ID and GOs
awk -F'\t' 'NR>1 {print $1, $10}' OFS='\t' ec_eggnog.tsv > eggnog_collapsed.tsv




awk -F'\t' '
{
    split($2, arr, ",")
    for (i in arr) {
        go[$1][arr[i]] = 1   # nested array: gene -> GO term -> exists
    }
}
END {
    for (g in go) {
        n = 0
        for (t in go[g]) {
            terms[n++] = t
        }
        asort(terms)
        out = terms[1]
        for (i=2; i<=n; i++) out = out "," terms[i]
        # strip .tN suffix from gene id here before printing
        gene = g
        sub(/\.t[0-9]+$/, "", gene)
        print gene "\t" out
        delete terms
    }
}
' blast_go_collapsed.tsv ipr_collapsed.tsv eggnog_collapsed.tsv | sort -k1,1 > merged_output.tsv


awk -F'\t' '{
    split($2, arr, ",")
    for (i in arr) {
        print $1 "\t" arr[i]
    }
}' merged_output.tsv > merged_output_long.tsv


