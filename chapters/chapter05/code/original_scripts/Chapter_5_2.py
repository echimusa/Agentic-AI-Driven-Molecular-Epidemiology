import io
import urllib.request
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# 1. Access aligned genomic data
data_url = "https://raw.githubusercontent.com/colinmeir/molecular_epi/main/data/influenza_h3n2_outbreak.fasta"
with urllib.request.urlopen(data_url) as response:
    fasta_data = response.read().decode('utf-8')

# Parse alignment with Biopython
alignment = AlignIO.read(io.StringIO(fasta_data), "fasta")
print(f"Alignment loaded successfully: {len(alignment)} samples, {alignment.get_alignment_length()} positions")

# 2. Calculate Distance Matrix using the Identity (JC-like) model
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)
print("\nFirst few lines of computed genetic distance matrix:")
print(str(distance_matrix)[:250] + " ...")

# 3. Construct a Neighbor-Joining Tree
constructor = DistanceTreeConstructor(calculator, 'nj')
nj_tree = constructor.build_tree(alignment)

# Label tree node supports based on standard nomenclature
print("\nTree construction finished. Tree Topology:")
Phylo.draw_ascii(nj_tree)

# Print root-to-tip genetic distances for individual samples
for leaf in nj_tree.get_terminals():
    print(f"Leaf {leaf.name:<15} distance to root: {nj_tree.distance(nj_tree.root, leaf):.5f}")