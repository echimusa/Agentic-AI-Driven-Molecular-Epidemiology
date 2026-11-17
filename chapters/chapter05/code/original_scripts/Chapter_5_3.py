import io
import os
import urllib.request
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from crewai import Agent, Task, Crew, Process

# 1. Define Agent Roles and Skills
# A specialized Agent that validates raw sequence alignments
dna_data_agent = Agent(
    role="Sequence Alignment & QC Specialist",
    goal="Ensure the multiple sequence alignment is high quality, has sufficient coverage, and is ready for inference.",
    backstory="You are an expert bioinformatician. You inspect multiple sequence alignments (MSAs), analyze GC content, identify and mask uncalled bases (N's), and guarantee that only robust sequence alignments proceed to evolutionary modeling.",
    verbose=True,
    allow_delegation=False
)

# A specialized Agent that evaluates evolutionary models
model_select_agent = Agent(
    role="Evolutionary Modeling & Statistical Officer",
    goal="Test and identify the mathematically optimal nucleotide substitution model (e.g., JC69, K80, or GTR) for the dataset.",
    backstory="You are a biomathematician. You evaluate nucleotide alignments to choose the best substitution model based on the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC), ensuring the analysis is biologically accurate.",
    verbose=True,
    allow_delegation=False
)

# An Agent that runs and interprets the tree
phylogenetics_agent = Agent(
    role="Phylogenetic Reconstruction Engine & Outbreak Investigator",
    goal="Reconstruct evolutionary trees, analyze transmission pathways, and identify outbreak clusters.",
    backstory="You are a molecular epidemiologist. You build phylogenetic trees (e.g., using Neighbor-Joining or Maximum Likelihood), assess branch confidence (bootstrap support), and trace transmission pathways to identify point-source outbreaks.",
    verbose=True,
    allow_delegation=False
)

# 2. Assign Sequential Tasks
# Task A: Validate and inspect the alignment
qc_task = Task(
    description=(
        "Review the sequence alignment. Download and inspect the FASTA file located at: "
        "'https://raw.githubusercontent.com/colinmeir/molecular_epi/main/data/influenza_h3n2_outbreak.fasta'. "
        "Summarize the sequence count, average read length, and ensure the data contains no invalid characters."
    ),
    expected_output="A quality control report confirming sequence count, alignment length, and dataset viability.",
    agent=dna_data_agent
)

# Task B: Select the evolutionary model
modeling_task = Task(
    description=(
        "Analyze the validated alignment. Determine whether a simple Jukes-Cantor (JC69) model is sufficient "
        "or if a complex, parameterized model like the Kimura 2-Parameter (K80) is required "
        "due to imbalances between transition and transversion rates."
    ),
    expected_output="A modeling report recommending the best substitution model with statistical justification.",
    agent=model_select_agent
)

# Task C: Build the tree and analyze transmission
tree_task = Task(
    description=(
        "Run the phylogenetic reconstruction. Use Biopython's DistanceCalculator with the recommended model "
        "to construct a Neighbor-Joining (NJ) tree from the influenza dataset. Interpret the resulting topology, "
        "and list any closely related pairs (genetic distance < 0.005) that suggest active transmission."
    ),
    expected_output="An epidemiological report describing the tree's structure, identifying clusters, and naming potential transmission pairs.",
    agent=phylogenetics_agent
)

# 3. Initialize and Execute the AI Flow
outbreak_analysis_crew = Crew(
    agents=[dna_data_agent, model_select_agent, phylogenetics_agent],
    tasks=[qc_task, modeling_task, tree_task],
    process=Process.sequential,  # Run tasks sequentially
    verbose=True
)

print("🚀 Starting the AI-Agentic Phylogenetic Analysis Workflow...")
result = outbreak_analysis_crew.kickoff()

print("\n================== FINAL AI-AGENT ANALYSIS OUTPUT ==================\n")
print(result)