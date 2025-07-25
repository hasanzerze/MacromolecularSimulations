#!/bin/bash

# Set fixed parameters
LATTICE_SIZE=500
EPSILON=-1.0
NUM_TRIALS=20000
SEED=42

# Output file
OUTPUT="ree_vs_n_eps-1.dat"
echo "# N  Ree2  Rg2" > $OUTPUT

# Sweep over chain lengths
for N in {10..80..10}; do
    echo "Running for CHAIN_LENGTH = $N"

    # Write input.dat
    cat << EOF > input.dat
LATTICE_SIZE $LATTICE_SIZE
CHAIN_LENGTH $N
EPSILON $EPSILON
NUM_TRIALS $NUM_TRIALS
SEED $SEED
EOF

    # Run the simulation
    ./rosenbluth > tmp_output.txt

    # Extract Ree^2 and Rg^2
    ree2=$(grep "Unbiased average Ree2" tmp_output.txt | awk '{print $NF}')
    rg2=$(grep "Unbiased average Rg2" tmp_output.txt | awk '{print $NF}')

    # Append to output
    echo "$N  $ree2  $rg2" >> $OUTPUT
done

echo "Done! Results saved to $OUTPUT"

