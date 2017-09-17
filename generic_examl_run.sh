#/usr/bin/env bash

# By Ryan Folk 2016, updated 2017

# Script for examl bootstrapping
q="partition.txt" # Partition
s="concatenated_protein_aware.phy" # Alignment
t=35 # Threads

# Executable paths
r="raxmlHPC-PTHREADS-SSE3"  # Old versions won't give specific partition files
e="examl-AVX"
p="parse-examl"


# VERY IMPORTANT
# Reorder by codon position using the individual position alignments generated during protein aware alignment
# To allow for stripping sites while maintaining reading frame.

# Make sure a new version of RAxML is installed that outputs bootstrap-specific partition files
# Make sure it is doing the right thing with codon features; reordering sites is normal behavior

# Create replicate datasets and partition files
# At this point any remaining undetermined sites will be removed
${r} -N 100 -b $RANDOM -f j -m GTRGAMMA -s ${s} -q ${q} -n rep

# Examl options: D cuts off likelihood search based on RF distance
for (( i=0 ; i<=99; i++ )); # Iterate through replicate datasets
do 
${p} -s ./${s}.BS${i} -m DNA -q ${q}.BS${i} -n parsed${i}
${r} -y -m GTRGAMMA -p $RANDOM -s ./${s}.BS${i} -n parsimony${i}
mpirun -np ${t} ${e} -t ./RAxML_parsimonyTree.parsimony${i} -m GAMMA -s ./parsed${i}.binary -D -n bootstrap${i}
done

# Concatenate all bootstrap replicates
cat ExaML_result.bootstrap* > ExaML_result.allbootstraps
cat RAxML_parsimonyTree.parsimony* > RAxML_parsimonyTree.allparsimony


# Calculate consensus of bootstraps for QC
${r} -J MR -m GTRGAMMA -s ${s} -z ExaML_result.allbootstraps -n majority_rule
${r} -J MR -m GTRGAMMA -s ${s} -z RAxML_parsimonyTree.allparsimony -n parsimony_check


# Remove useless branch lengths and non-standard newick format

sed -i 's/:1\.0\[//g' RAxML_MajorityRuleConsensusTree.majority_rule
sed -i 's/\]//g' RAxML_MajorityRuleConsensusTree.majority_rule
sed -i 's/:1\.0\[//g' RAxML_MajorityRuleConsensusTree.parsimony_check
sed -i 's/\]//g' RAxML_MajorityRuleConsensusTree.parsimony_check

# Examl best tree search, use every 10th bootstrap as a starting tree like so:
# The RF cutoff is not used here
# Iteration is correct
${p} -s ${s} -m DNA -q ${q} -n parsed_besttree
for (( i=0 ; i<=99; i+=10 )); # Iterate through replicate datasets
do 
mpirun -np ${t} ${e} -t ./RAxML_parsimonyTree.parsimony${i} -m GAMMA -s ./parsed_besttree.binary -n besttree${i}
done

grep "Likelihood   : " ExaML_info.besttree* | sed 's/besttree0/besttree00/g' # Take maximum

# Plot on best tree
z=    #INSERT BEST NUMBER HERE based on reported likelihood by grep
${r} -f b -m GTRGAMMA -s ${s} -z ExaML_result.allbootstraps -t ExaML_result.besttree${z} -n bootstrap_besttree

rm ExaML_binaryCheckpoint*
rm *binary

for (( i=0 ; i<=99; i++ )); # Iterate through replicate datasets
do
rm *.BS${i}
done
