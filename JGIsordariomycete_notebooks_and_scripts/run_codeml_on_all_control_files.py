import os
import subprocess
import shlex

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

def run_codeml(control_file_path):
    print(f"Running codeml on: {control_file_path}")

    # Check if the control file exists
    if not os.path.exists(control_file_path):
        print(f"Control file not found: {control_file_path}")
        return
    
    # Verify and print the contents of the control file
    with open(control_file_path, 'r') as file:
        control_file_content = file.read()
        print(f"Control file contents:\n{control_file_content}")

    # Explicitly verify the seqfile and treefile paths in the control file
    control_lines = control_file_content.splitlines()
    seqfile_path = None
    treefile_path = None

    for line in control_lines:
        if line.startswith("seqfile"):
            seqfile_path = line.split("=")[1].strip()
            print(f"Seqfile path: {seqfile_path}")
        if line.startswith("treefile"):
            treefile_path = line.split("=")[1].strip()
            print(f"Treefile path: {treefile_path}")

    if seqfile_path and not os.path.exists(seqfile_path):
        print(f"Sequence file not found: {seqfile_path}")
        return
    if treefile_path and not os.path.exists(treefile_path):
        print(f"Tree file not found: {treefile_path}")
        return

    # Run codeml directly with detailed logging
    log_file = control_file_path.replace(".ctl", ".log")
    command = f"codeml {control_file_path}"
    print(f"Running command: {command}")
    with open(log_file, 'w') as log:
        process = subprocess.Popen(shlex.split(command), stdout=log, stderr=log)
        try:
            stdout, stderr = process.communicate(timeout=300)  # Timeout after 300 seconds
            print(f"stdout: {stdout.decode() if stdout else 'None'}")
            print(f"stderr: {stderr.decode() if stderr else 'None'}")
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
            print(f"Process timed out. Killed process for control file: {control_file_path}")
            print(f"stdout: {stdout.decode() if stdout else 'None'}")
            print(f"stderr: {stderr.decode() if stderr else 'None'}")
        
    print(f"Successfully ran codeml on {control_file_path}")

def process_files(control_dir):
    for filename in os.listdir(control_dir):
        if filename.endswith(".ctl"):
            control_file_path = os.path.join(control_dir, filename)
            print(f"\nProcessing control file: {control_file_path}")
            run_codeml(control_file_path)
            break  # Process one file at a time for detailed debugging

# Example usage
control_directory = '/usr2/people/shollyt22/shollyt22/JGIsordariomycete/after_conversion/control_files_nul'
process_files(control_directory)
