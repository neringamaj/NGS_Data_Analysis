import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
import subprocess
import pandas as pd

fastq_file_path = 'reads_for_analysis.fastq'

def read_fastq_scores(fastq_file_path):
    quality_scores = []
    sequences = []
    with open(fastq_file_path, 'r') as file:
        while True:
            identifier = file.readline().strip() 
            if not identifier:
                break
            sequence = file.readline().strip() 
            file.readline()  # Skip the '+'
            quality_scores.append(file.readline().strip()) 
            sequences.append((identifier, sequence))
    return quality_scores, sequences

def determine_encoding(quality_scores):
    encoding_ranges = {
        'Sanger Phred+33': range(32, 78),
        'Solexa Solexa+64': range(58, 105),
        'Illumina 1.3+ Phred+64': range(63, 105),
        'Illumina 1.5+ Phred+64': range(65, 106),
        'Illumina 1.8+ Phred+33': range(32, 75)
    }
    
    possible_encodings = {encoding: True for encoding in encoding_ranges}
    
    for line in quality_scores:
        for char in line:
            ascii_value = ord(char)
            for encoding, ascii_range in encoding_ranges.items():
                if ascii_value not in ascii_range:
                    possible_encodings[encoding] = False
    
    return [encoding for encoding, possible in possible_encodings.items() if possible]

def calculate_cg_proportion(sequence):
    cg_count = sequence.count('C') + sequence.count('G')
    return cg_count / len(sequence) if sequence else 0

def run_blast(query_file, database, output_file, max_seqs, entrez_query):

    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames"

    blast_command = [
        'blastn' if database == 'nt' else 'blastx', 
        '-query', query_file, 
        '-db', database, 
        '-remote',
        '-out', output_file,
        '-outfmt', outfmt,
        '-entrez_query', entrez_query, 
        '-max_target_seqs', str(max_seqs),
    ]
    try:
        subprocess.run(blast_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running BLAST: {e}")

def read_blast_output_and_create_table(output_file_nt, output_file_nr):
    nt_results = pd.read_csv(output_file_nt, header=None, delimiter='\t')
    nr_results = pd.read_csv(output_file_nr, header=None, delimiter='\t')
    
    nt_results = nt_results[[0, 12, 13]]
    nr_results = nr_results[[0, 12, 13]]
    
    nt_results.columns = ['Seq ID', 'Read ID', 'Species']
    nr_results.columns = ['Seq ID', 'Read ID', 'Species']
    
    combined_results = pd.concat([nt_results, nr_results])

    return combined_results

def write_sequences_to_fasta(sequences, filepath):
    with open(filepath, 'w') as file:
        for seq_id, seq in sequences:
            seq_id_formatted = seq_id.split()[0].replace(' ', '_')
            file.write(f">{seq_id_formatted}\n{seq}\n")

def plot_table(table_data):
    row_height = 0.5
    fig_height = len(table_data) * row_height
    fig, ax = plt.subplots(figsize=(12, fig_height))
    ax.axis('tight')
    ax.axis('off')
    table = ax.table(cellText=table_data.values,
                     colLabels=table_data.columns,
                     colWidths=[0.6,0.15,0.35],
                     loc='center',
                     cellLoc='center')
    cell_height = table.get_celld()[0,0].get_height()
    for key, cell in table.get_celld().items():
        cell.set_height(cell_height * 2)
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 1)  
    plt.savefig('table.png')
    plt.show()

# Main logic
fastq_file_path = 'reads_for_analysis.fastq'

quality_score_lines, sequences = read_fastq_scores(fastq_file_path)
possible_encodings = determine_encoding(quality_score_lines)
print(f"Possible quality encodings: {possible_encodings}")

identifiers, seqs = zip(*sequences) 
cg_proportions = np.array([calculate_cg_proportion(seq) for seq in seqs])

num_bins = 50  
hist_counts, bin_edges = np.histogram(cg_proportions, bins=num_bins, range=(0,1))
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

plt.figure(figsize=(10, 6))
plt.bar(bin_centers, hist_counts, width=1/num_bins, color='skyblue')
plt.title('Distribution of C/G Nucleotide Proportions in Read Sequences')
plt.xlabel('Proportion of C/G Nucleotides')
plt.ylabel('Number of Reads')
plt.grid(True)
plt.savefig('cg.png')
plt.show()

peak_threshold = 0.75 * np.max(hist_counts)
peaks, _ = find_peaks(hist_counts, height=peak_threshold)

selected_sequences = {}
for peak in peaks:
    peak_start = bin_edges[peak]
    peak_end = bin_edges[peak + 1]
    sequences_in_peak = [(identifiers[i], seqs[i]) for i, cg in enumerate(cg_proportions) if peak_start <= cg < peak_end]
    top_sequences = sorted(sequences_in_peak, key=lambda x: x[1], reverse=True)[:5]
    selected_sequences[peak] = top_sequences

all_selected_sequences = [seq for sequences in selected_sequences.values() for seq in sequences]

output_fasta_path = 'selected_sequences.fasta'
write_sequences_to_fasta(all_selected_sequences, output_fasta_path)

print(f"Selected sequences written to {output_fasta_path}")

run_blast('selected_sequences.fasta', 'nt', 'blast_output_nt.csv', 1, 'bacteria[organism]')
run_blast('selected_sequences.fasta', 'nr', 'blast_output_nr.csv', 1, 'bacteria[organism]')

table_csv = read_blast_output_and_create_table('blast_output_nt.csv', 'blast_output_nr.csv')
plot_table(table_csv)