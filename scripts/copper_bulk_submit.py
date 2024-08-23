import os
import subprocess
from pathlib import Path, PosixPath


def grant_execute_permissions(directory_path) -> None:
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


def _get_ssh_username() -> str:
    try:
        return os.environ["HPC_USERNAME"]
    except KeyError:
        return input("Username: ")


def _get_hpc_workspace_directory() -> PosixPath:
    try:
        return PosixPath(os.environ["HPC_WORKSPACE"])
    except KeyError:
        username = _get_ssh_username()
        return PosixPath(f"/rds/user/{username}/hpc-work")


def _get_relative_hpc_path(remote_folder: PosixPath) -> PosixPath:
    workspace_folder = _get_hpc_workspace_directory()
    return workspace_folder / remote_folder


def copy_files_to_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files from local to remote."""
    username = _get_ssh_username()

    remote_folder_absolute = _get_relative_hpc_path(remote_folder)
    command = [
        "scp",
        "-r",
        f"{local_folder.absolute()}",
        f"{username}@login.hpc.cam.ac.uk:{remote_folder_absolute.as_posix()}",
    ]
    subprocess.run(
        command,  # noqa: S603
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        check=True,
    )


def copy_files_from_hpc(local_folder: Path, remote_folder: PosixPath) -> None:
    """Make use of the scp command to copy files from local to remote."""
    username = _get_ssh_username()

    remote_folder_absolute = _get_relative_hpc_path(remote_folder)
    command = [
        "scp",
        "-r",
        f"{username}@login.hpc.cam.ac.uk:{remote_folder_absolute.as_posix()}",
        f"{local_folder.absolute()}",
    ]
    subprocess.run(
        command,  # noqa: S603
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        check=True,
    )


def run_command_on_hcp(remote_command: str) -> None:
    """Use SSH to run remote_command on the hpc login node."""
    username = _get_ssh_username()
    # Works in terminal, but not in here
    command = [
        "ssh",
        "-f",
        f"{username}@login.hpc.cam.ac.uk",
        f"<<ENDSSH\n{remote_command}\nENDSSH",
    ]
    subprocess.run(
        command,  # noqa: S603
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        check=True,
    )


if __name__ == "__main__":
    # Define the paths
    local_path = "./data"
    remote_path = "ll720@login.hpc.cam.ac.uk:/rds/user/ll720/hpc-work"

    # copy_files_to_hpc(Path("./data"), PosixPath("./data"))
    command = """
cd /rds/user/mo433/hpc-work
touch test.test
"""
    run_command_on_hcp(command)

    # # Grant execute permissions to all files in the local directory
    # grant_execute_permissions(local_path)

    # # Run the scp command
    # run_scp_command(local_path, remote_path)
