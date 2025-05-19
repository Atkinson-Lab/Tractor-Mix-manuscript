import sys
from collections import defaultdict

def find_duplicate_variant_ids(bim_file, output_file):
    # Dictionary to store variant IDs and their line numbers
    variant_dict = defaultdict(list)
    
    # Reading the .bim file and storing variant IDs
    with open(bim_file, 'r') as f:
        for line_number, line in enumerate(f, 1):
            columns = line.strip().split()
            variant_id = columns[1]  # Variant ID is the second column in the .bim file
            variant_dict[variant_id].append(line_number)
    
    # Writing duplicates to the output file
    with open(output_file, 'w') as f_out:
        for variant_id, lines in variant_dict.items():
            if len(lines) > 1:
                f_out.write(f"{variant_id}\n")
    
    print(f"Duplicate variant IDs written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python find_duplicates.py <input.bim> <output.txt>")
        sys.exit(1)
    
    bim_file = sys.argv[1]
    output_file = sys.argv[2]
    
    find_duplicate_variant_ids(bim_file, output_file)
