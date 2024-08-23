import os
import subprocess


def grant_execute_permissions(directory) -> None:
    """Grant execute permissions to all files in the given directory."""
    for root, _dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            os.chmod(file_path, 0o755)  # 0o755 gives execute permissions


def run_scp_command(local_path, remote_path, username) -> None:
    """Run the SCP command to copy files from local to remote."""
    scp_command = f"scp -r {local_path} {remote_path}"
    try:
        subprocess.run(scp_command, shell=True, check=True)
        print(f"Files successfully copied to {remote_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error during SCP transfer: {e}")


def submit_hpc_job():
    """Submit the HPC job using an appropriate command and return the job ID."""
    # Example HPC submission command; replace with actual command for your HPC system
    submit_command = "sbatch submit_job_script.sh"
    try:
        result = subprocess.run(
            submit_command, shell=True, check=True, capture_output=True
        )
        print(f"Job submitted: {result.stdout.decode().strip()}")
        return (
            result.stdout.decode().strip().split()[-1]
        )  # Extract job ID from the output
    except subprocess.CalledProcessError as e:
        print(f"Error submitting HPC job: {e}")
        return None


def load_hpc_result(job_id) -> None:
    """Load and process results from the HPC job."""
    # Example: Check the status of the job and retrieve results
    if job_id:
        print(f"Loading results for job ID: {job_id}")
        # Add logic to load results
    else:
        print("No job ID provided; cannot load results.")
