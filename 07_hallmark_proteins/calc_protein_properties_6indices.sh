#!/bin/bash

# check the number of parameters
if [[ $# -ne 2 ]]; then
    echo "usage: $0 <input faa file> <output file>"
    exit 1
fi

# input and output
input_file="$1"
output_file="$2"

# setting amino acid scales
declare -A hydropathicity=(
    [A]=1.8   [R]=-4.5  [N]=-3.5  [D]=-3.5  [C]=2.5
    [Q]=-3.5  [E]=-3.5  [G]=-0.4  [H]=-3.2  [I]=4.5
    [L]=3.8   [K]=-3.9  [M]=1.9   [F]=2.8   [P]=-1.6
    [S]=-0.8  [T]=-0.7  [W]=-0.9  [Y]=-1.3  [V]=4.2
)

declare -A relative_mutability=(
    [A]=100.0 [R]=65.0  [N]=134.0 [D]=106.0 [C]=20.0
    [Q]=93.0  [E]=102.0 [G]=49.0  [H]=66.0  [I]=96.0
    [L]=40.0  [K]=56.0  [M]=94.0  [F]=41.0  [P]=56.0
    [S]=120.0 [T]=97.0  [W]=18.0  [Y]=41.0  [V]=74.0
)

declare -A average_flexibility=(
    [A]=0.360 [R]=0.530 [N]=0.460 [D]=0.510 [C]=0.350
    [Q]=0.490 [E]=0.500 [G]=0.540 [H]=0.320 [I]=0.460
    [L]=0.370 [K]=0.470 [M]=0.300 [F]=0.310 [P]=0.510
    [S]=0.510 [T]=0.440 [W]=0.310 [Y]=0.420 [V]=0.390
)

declare -A refractivity=(
    [A]=4.34  [R]=26.66 [N]=13.28 [D]=12.00 [C]=35.77
    [Q]=17.56 [E]=17.26 [G]=0.00  [H]=21.81 [I]=19.06
    [L]=18.78 [K]=21.29 [M]=21.64 [F]=29.40 [P]=10.93
    [S]=6.35  [T]=11.01 [W]=42.53 [Y]=31.53 [V]=13.92
)

declare -A polarity=(
    [A]=8.10  [R]=10.50 [N]=11.60 [D]=13.00 [C]=5.50
    [Q]=10.50 [E]=12.30 [G]=9.00  [H]=10.40 [I]=5.20
    [L]=4.90  [K]=11.30 [M]=5.70  [F]=5.20  [P]=8.00
    [S]=9.20  [T]=8.60  [W]=5.40  [Y]=6.20  [V]=5.90
)

declare -A transmembrane_tendency=(
    [A]=0.38  [R]=-2.57 [N]=-1.62 [D]=-3.27 [C]=-0.30
    [Q]=-1.84 [E]=-2.90 [G]=-0.19 [H]=-1.44 [I]=1.97
    [L]=1.82  [K]=-3.46 [M]=1.40  [F]=1.98  [P]=-1.44
    [S]=-0.53 [T]=-0.32 [W]=1.53  [Y]=0.49  [V]=1.46
)

# clear output file and write the header
echo -e "Sequence_ID\tLength\tHydropathicity\tRelative_mutability\tAverage_flexibility\tRefractivity\tPolarity\tTransmembrane_tendency" > "$output_file"

# use awk to calculate 6 indices of each protein
awk -v RS='>' 'NR > 1 {
    split($0, lines, "\n")
    id = lines[1]
    seq = ""
    for (i = 2; i <= length(lines); i++) {
        gsub("\\*", "", lines[i])  #remove asterisk
        seq = seq lines[i]
    }
    print id "\t" seq
}' "$input_file" | while IFS=$'\t' read -r id seq; do
    
    len=${#seq}

    if [ "$len" -eq 0 ]; then
        echo -e "$id\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "$output_file"
        continue
    fi

    gravy_total=0
    polarity_total=0
    relative_mutability_total=0
    average_flexibility_total=0
    refractivity_total=0
    transmembrane_tendency_total=0

    for ((i=0; i<len; i++)); do
        aa=${seq:i:1}

        # hydropathicity
        g=${hydropathicity[$aa]}
        [[ -z $g ]] && g=0
        gravy_total=$(echo "$gravy_total + $g" | bc)

        # relative_mutability
        relm=${relative_mutability[$aa]}
        [[ -z $relm ]] && relm=0
        relative_mutability_total=$(echo "$relative_mutability_total + $relm" | bc)
		
	    # average_flexibility
        af=${average_flexibility[$aa]}
        [[ -z $af ]] && af=0
        average_flexibility_total=$(echo "$average_flexibility_total + $af" | bc)
		
	    # refractivity
        ref=${refractivity[$aa]}
        [[ -z $ref ]] && ref=0
        refractivity_total=$(echo "$refractivity_total + $ref" | bc)
		
	    # polarity
        p=${polarity[$aa]}
        [[ -z $p ]] && p=0
        polarity_total=$(echo "$polarity_total + $p" | bc)
		
	    # transmembrane_tendency
        tt=${transmembrane_tendency[$aa]}
        [[ -z $tt ]] && tt=0
        transmembrane_tendency_total=$(echo "$transmembrane_tendency_total + $tt" | bc)
    done

    gravy_avg=$(printf "%.3f" "$(echo "scale=5; $gravy_total / $len" | bc)")
    relative_mutability_avg=$(printf "%.3f" "$(echo "scale=5; $relative_mutability_total / $len" | bc)")
    average_flexibility_avg=$(printf "%.3f" "$(echo "scale=5; $average_flexibility_total / $len" | bc)")
    refractivity_avg=$(printf "%.3f" "$(echo "scale=5; $refractivity_total / $len" | bc)")
    polarity_avg=$(printf "%.3f" "$(echo "scale=5; $polarity_total / $len" | bc)")
    transmembrane_tendency_avg=$(printf "%.3f" "$(echo "scale=5; $transmembrane_tendency_total / $len" | bc)")

    echo -e "${id}\t${len}\t${gravy_avg}\t${relative_mutability_avg}\t${average_flexibility_avg}\t${refractivity_avg}\t${polarity_avg}\t${transmembrane_tendency_avg}" >> "$output_file"
done

