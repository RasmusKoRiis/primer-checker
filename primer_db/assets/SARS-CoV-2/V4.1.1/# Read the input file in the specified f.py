# Read the input file in the specified format
input_file_path = "/home/rasmuskopperud.riis/Coding/run_folder/nf-core-sarscovseq/assets/V4.1/SARS-CoV-2.primer.bed"  # Replace with your actual input file path
output_file_path = "/home/rasmuskopperud.riis/Coding/run_folder/nf-core-sarscovseq/assets/V4.1/primers.fasta"

with open(input_file_path, "r") as infile, open(output_file_path, "w") as outfile:
    for line in infile:
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue  # Skip any malformed lines
        identifier = parts[3]   # Assuming the fourth column is the identifier
        sequence = parts[6]     # Assuming the seventh column is the sequence

        # Write in FASTA format
        outfile.write(f">{identifier}\n{sequence}\n")

print(f"Parsing complete. Output saved to {output_file_path}")
