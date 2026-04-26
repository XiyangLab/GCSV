#!/bin/bash

# ============================================================
# Step 01: Sequence similarity search against PhaLP database
# Tool: MMseqs2(v18-8cc5c)
# Description:
#   - Build databases for reference and query sequences
#   - Perform sequence similarity search
#   - Extract best hit (lowest e-value) per query
#
# Input:
#   - QUERY_FASTA: protein sequences (FASTA)
#   - REFERENCE_FASTA: proteins from PhaLP database
#
# Output:
#   - final_results.m8: raw alignment results
#   - best_hits.csv: best hit per query (filtered)
# ============================================================

# ========== Configuration ==========
THREADS=16
EVALUE=1e-5

QUERY_FASTA="input/query.faa"
REFERENCE_FASTA="input/reference_phalp.faa"

WORKDIR="mmseqs_work"
TMPDIR="${WORKDIR}/tmp"

QUERY_DB="${WORKDIR}/queryDB"
REF_DB="${WORKDIR}/refDB"
RESULT_DB="${WORKDIR}/resultDB"

ALIGN_RESULT="${WORKDIR}/final_results.m8"
BEST_HITS="${WORKDIR}/best_hits.csv"

mkdir -p ${WORKDIR} ${TMPDIR}

# ========== Step 1: Create databases ==========
echo "Creating MMseqs2 databases..."

mmseqs createdb ${REFERENCE_FASTA} ${REF_DB}
mmseqs createdb ${QUERY_FASTA} ${QUERY_DB}

# ========== Step 2: Sequence search ==========
echo "Running MMseqs2 search..."

mmseqs search ${QUERY_DB} \
              ${REF_DB} \
              ${RESULT_DB} \
              ${TMPDIR} \
              -e ${EVALUE} \
              --threads ${THREADS} \
              -s 5.5

# ========== Step 3: Convert results ==========
echo "Converting alignment results..."

mmseqs convertalis ${QUERY_DB} \
                   ${REF_DB} \
                   ${RESULT_DB} \
                   ${ALIGN_RESULT}

# ========== Step 4: Extract best hits ==========
echo "Filtering best hits (lowest e-value per query)..."

awk '
function sci_to_num(s) {
    if (match(s, /([0-9.]+)[Ee]([+-]?[0-9]+)/, arr)) {
        return arr[1] * 10^arr[2]
    }
    return s
}
{
    query_id = $1
    evalue_str = $11
    evalue_num = sci_to_num(evalue_str)

    if (!(query_id in min_evalue) || evalue_num < min_evalue[query_id]) {
        min_evalue[query_id] = evalue_num
        best_line[query_id] = $0
    }
}
END {
    print "query_id,target_id,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"
    for (qid in best_line) {
        print best_line[qid] | "tr \"\t\" \",\""
    }
}
' ${ALIGN_RESULT} > ${BEST_HITS}

echo "Done!"
echo "Best hits saved to: ${BEST_HITS}"