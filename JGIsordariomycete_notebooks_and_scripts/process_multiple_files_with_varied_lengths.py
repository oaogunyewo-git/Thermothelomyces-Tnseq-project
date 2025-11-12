import os
import re

def clean_and_correct_sequence_data(file_content, expected_length):
    corrected_content = []
    for line in file_content:
        if re.match(r'^[a-zA-Z]', line):
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                name, sequence = parts
                # Remove any invalid characters (anything not A, T, G, C, N, or -)
                sequence = re.sub(r'[^ATGCN-]', '', sequence.upper())
                # Pad the sequence if it's shorter than the expected length
                sequence = sequence.ljust(expected_length, '-')
                corrected_content.append(f"{name}  {sequence}\n")
            else:
                corrected_content.append(line)
        else:
            corrected_content.append(line)
    return corrected_content

def verify_and_fix_sequence_file(file_path):
    print(f"Processing file: {file_path}")
    with open(file_path, 'r') as input_file:
        file_content = input_file.readlines()
    
    # Determine the expected length by finding the longest sequence
    expected_length = 0
    for line in file_content:
        if re.match(r'^[a-zA-Z]', line):
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                sequence = parts[1]
                expected_length = max(expected_length, len(sequence))
    
    print(f"Expected length for sequences in {file_path}: {expected_length}")
    
    corrected_content = clean_and_correct_sequence_data(file_content, expected_length)
    
    with open(file_path, 'w') as corrected_file:
        corrected_file.writelines(corrected_content)
    
    print(f"Corrected file saved to: {file_path}")
    return file_path

def create_control_file(seqfile, treefile, outfile, control_file_path):
    control_content = f"""
seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 9
verbose = 0
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0
model = 2

NSsites = 0
icode = 0
Mgene = 0

fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 2

fix_alpha = 1
alpha = .0
Malpha = 0
ncatG = 3

getSE = 0
RateAncestor = 0

fix_blength = 0
method = 0
Small_Diff = .45e-6
cleandata = 1
"""

    # Remove leading and trailing spaces from each line
    control_content = "\n".join(line.strip() for line in control_content.strip().split("\n"))

    with open(control_file_path, 'w') as control_file:
        control_file.write(control_content)

    print(f"Control file created: {control_file_path}")

def process_files(input_dir, output_dir, tree_dir, control_dir):
    # Ensure output directories exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(control_dir):
        os.makedirs(control_dir)

    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".phylip"):
            input_file_path = os.path.join(input_dir, filename)

            # Skip directories
            if os.path.isdir(input_file_path):
                continue

            try:
                # Read and correct the sequence file
                corrected_file_path = verify_and_fix_sequence_file(input_file_path)

                # Determine tree and output file paths
                treefile = os.path.join(tree_dir, filename.replace(".phylip", ".phylip_pruned.nwk"))
                if not os.path.exists(treefile):
                    print(f"Tree file not found: {treefile}")
                    continue

                outfile = os.path.join(output_dir, filename.replace(".phylip", ".altout"))

                # Create control file
                control_file_path = os.path.join(control_dir, filename.replace(".phylip", ".ctl"))
                create_control_file(corrected_file_path, treefile, outfile, control_file_path)

                print(f"Processed {filename}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

# Example usage
input_directory = '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/after_conversion'
output_directory = '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/For_second_batch_singlecopy_orthogroup_May_09_2024/written_OG_nucl/Filtered_OG_divisible_by_3/paml_output_files_alt/from_MACSE_files_after_conversion'
tree_directory = '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/For_second_batch_singlecopy_orthogroup_May_09_2024/written_OG_nucl/Filtered_OG_divisible_by_3/Pruned_newick_trees_alt'
control_directory = '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/after_conversion/control_files'

process_files(input_directory, output_directory, tree_directory, control_directory)
