# The final report generated autonomously by the agents:
# "Dataset contains 15 aligned Influenza H3N2 sequences (975 bp). No corrupt bases were detected.
# Jukes-Cantor (JC69) was selected based on genetic uniformity across early sampling periods.
# Top pairs identified: Strain_A and Strain_B (Genetic Distance = 0.0012) and Strain_G and Strain_H (0.0031). 
# We recommend immediate epidemiological follow-up for Patients A, B, G, and H to investigate shared exposures."

# Verify by running the model and visualizing the tree:
with urllib.request.urlopen("https://raw.githubusercontent.com/colinmeir/molecular_epi/main/data/influenza_h3n2_outbreak.fasta") as res:
    alignment = AlignIO.read(io.StringIO(res.read().decode('utf-8')), "fasta")

calculator = DistanceCalculator('jc69')
dist_matrix = calculator.get_distance(alignment)
tree = DistanceTreeConstructor(calculator, 'nj').build_tree(alignment)

# Show final validation check
print(f"Direct Verification - Leaf nodes found: {len(tree.get_terminals())}")
for term in list(tree.get_terminals())[:3]:
    print(f" Verified Leaf ID: {term.name:<25} Root distance: {tree.distance(tree.root, term):.6f}")