import subprocess


def grant_execute_permissions(directory_path: str) -> None:
    """Grant execute permissions to all files in the given directory."""
    # Construct the chmod command
    command = ["chmod", "-R", "+x", directory_path]

    try:
        # Run the command
        result = subprocess.run(command, check=True, text=True, capture_output=True)

        # Print the command output
        print("Permission command output:", result.stdout)
    except subprocess.CalledProcessError as e:
        # Print the error message
        print("Error occurred while setting permissions:", e.stderr)


def run_scp_command(local_path, remote_path) -> None:
    """Run the scp command to copy files from local to remote."""
    # Construct the scp command
    command = ["scp", "-r", local_path, remote_path]

    try:
        # Run the command
        result = subprocess.run(command, check=True, text=True, capture_output=True)

        # Print the command output
        print("SCP command output:", result.stdout)
    except subprocess.CalledProcessError as e:
        # Print the error message
        print("Error occurred while transferring files:", e.stderr)
