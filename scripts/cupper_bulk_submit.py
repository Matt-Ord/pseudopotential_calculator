from pseudopotential_calculator import submit

# Define the paths
local_path = "./data"
remote_path = "ll720@login.hpc.cam.ac.uk:/rds/user/ll720/hpc-work"

if __name__ == "__main__":
    # Grant execute permissions to all files in the local directory
    submit.grant_execute_permissions(local_path)

    # Run the scp command
    submit.run_scp_command(local_path, remote_path)
