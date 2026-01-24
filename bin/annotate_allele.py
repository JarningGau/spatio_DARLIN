import pandas as pd
from darlinpy import analyze_sequences
import sys

locus = sys.argv[1] # CA,TA,RA
min_bc_len = sys.argv[2] 
in_feat = sys.argv[3]
out_feat_allele = sys.argv[4]
out_feat_original = sys.argv[5]

locus_dict = {
    "CA": "Col1a1",
    "TA": "Tigre",
    "RA": "Rosa"
}
locus_name = locus_dict[locus]

data = pd.read_csv(in_feat, sep="\t", header=None)
data.columns = ['query', 'clone_bc', 'type']
sequences = data['clone_bc'].tolist()
results_allele = analyze_sequences(
        sequences = sequences, 
        config = locus_name,
        min_sequence_length = int(min_bc_len), 
        verbose=False)
results_allele = results_allele.to_df()

data = pd.merge(data, results_allele[['query', 'mutations']], how="left", on="query")
data = data[['clone_bc', 'mutations', 'type']]
data.columns = ['clone_id', 'allele', 'type']

results_allele.to_csv(out_feat_original, sep="\t", index=False)
data.to_csv(out_feat_allele, sep="\t", index=False)
