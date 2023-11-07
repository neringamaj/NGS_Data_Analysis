import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

# Define the path to your FASTQ file here
fastq_file_path = 'reads_for_analysis.fastq'

# Function to read the quality score lines from the FASTQ file
def read_fastq_scores(fastq_file_path):
    quality_scores = []
    with open(fastq_file_path, 'r') as file:
        while True:
            identifier = file.readline().strip()  # Read the sequence identifier
            if not identifier:
                break
            file.readline()  # Skip the sequence itself
            file.readline()  # Skip the '+'
            quality_scores.append(file.readline().strip())  # Get the quality score line
    return quality_scores

# Function to determine the encoding from the ASCII values
def determine_encoding(quality_scores):
    encoding_ranges = {
        'Sanger Phred+33': range(33, 78),
        'Solexa Solexa+64': range(59, 105),
        'Illumina 1.3+ Phred+64': range(64, 105),
        'Illumina 1.5+ Phred+64': range(67, 105),
        'Illumina 1.8+ Phred+33': range(33, 75)
    }
    
    # Initialize the dictionary to hold possible encodings
    possible_encodings = {encoding: True for encoding in encoding_ranges}
    
    # Check each quality score line
    for line in quality_scores:
        for char in line:
            ascii_value = ord(char)
            # Update the possible encodings based on the ASCII value
            for encoding, ascii_range in encoding_ranges.items():
                if ascii_value not in ascii_range:
                    possible_encodings[encoding] = False
    
    # Return the encodings that still have True as their value
    return [encoding for encoding, possible in possible_encodings.items() if possible]

# Function to read the sequence lines from the FASTQ file
def read_fastq_sequences(fastq_file_path):
    sequences = []
    with open(fastq_file_path, 'r') as file:
        while True:
            identifier = file.readline().strip()  # Read the sequence identifier
            if not identifier:
                break
            sequence = file.readline().strip()  # Read the sequence itself
            file.readline()  # Skip the '+'
            file.readline()  # Skip the quality score line
            sequences.append((identifier, sequence))
    return sequences

# Function to calculate the CG proportion for a given sequence
def calculate_cg_proportion(sequence):
    # This function assumes 'sequence' is a string of nucleotides.
    cg_count = sequence.count('C') + sequence.count('G')
    return cg_count / len(sequence) if sequence else 0

# Function to read the sequence lines from the FASTQ file
def read_fastq_sequences(fastq_file_path):
    sequences = []
    with open(fastq_file_path, 'r') as file:
        while True:
            identifier = file.readline().strip()  # Read the sequence identifier
            if not identifier:
                break
            sequence = file.readline().strip()  # Read the sequence itself
            file.readline()  # Skip the '+'
            file.readline()  # Skip the quality score line
            sequences.append((identifier, sequence))
    return sequences

# Function to select sequences near peaks
def select_sequences_near_peaks(sequences, peaks, cg_proportions):
    selected_sequences = []
    for peak_idx in peaks:
        peak_cg_prop = cg_proportions[peak_idx]
        distances = np.abs(cg_proportions - peak_cg_prop)
        nearest_indices = np.argsort(distances)[:5]
        selected_sequences.extend([sequences[idx] for idx in nearest_indices])
    return selected_sequences

# Function to analyze CG content and find peaks
def analyze_cg_content_and_find_peaks(cg_proportions):
    peaks, _ = find_peaks(cg_proportions, height=0.85, distance=150)  # Adjust these parameters as needed

    print(f"Number of peaks found: {len(peaks)}")

    # If no peaks are found, adjust parameters or check data
    if len(peaks) == 0:
        print("No peaks found. Adjusting parameters for peak detection.")
        peaks, _ = find_peaks(cg_proportions, height=None, distance=1)  # Use default height and minimum distance

    # If more than 5 peaks are found, take only the 5 largest
    if len(peaks) > 5:
        peak_heights = cg_proportions[peaks]
        largest_peaks = np.argsort(peak_heights)[-5:]
        peaks = peaks[largest_peaks]

    # Print the number of peaks found
    print(cg_proportions[peaks])
    return peaks

# Main logic
fastq_file_path = 'reads_for_analysis.fastq'

# Check the quality encoding
quality_score_lines = read_fastq_scores(fastq_file_path)
possible_encodings = determine_encoding(quality_score_lines)
print(f"Possible quality encodings: {possible_encodings}")

sequences = read_fastq_sequences(fastq_file_path)
identifiers, seqs = zip(*sequences)  # Unzip the identifiers and sequences
cg_proportions = np.array([calculate_cg_proportion(seq) for seq in seqs])

# Find peaks and select sequences near those peaks
peaks = analyze_cg_content_and_find_peaks(cg_proportions)
selected_sequences = select_sequences_near_peaks(sequences, peaks, cg_proportions)

# Plotting the histogram of CG proportions
plt.figure(figsize=(10, 4))
plt.plot(cg_proportions, label='CG Proportion')
plt.scatter(peaks, cg_proportions[peaks], color='r', label='Peaks')
plt.title('C/G Nucleotide Proportion Distribution')
plt.xlabel('Sequence index')
plt.ylabel('CG Proportion')
plt.legend()
plt.savefig('cg_proportion.png')  # Save the histogram to a file
plt.close()  # Close the plot to free up memory

# Function to write sequences to a FASTA file
def write_sequences_to_fasta(sequences, filepath):
    with open(filepath, 'w') as file:
        for seq_id, seq in sequences:
            file.write(f">{seq_id}\n{seq}\n")

# Write selected sequences to a FASTA file
output_fasta_path = 'selected_sequences.fasta'
write_sequences_to_fasta(selected_sequences, output_fasta_path)

print(f"Selected sequences written to {output_fasta_path}")
