import os
import subprocess

import config  # Assuming config.py contains the range of n values


def run_castep_calculations(k_range) -> None:
    base_dir = "/workspaces/pseudopotential_calculator/data/"
    command_template = "castep.serial bulk_cu_{0}x{0}x{0}"

    for n in k_range:
        # Format the directory name and command with the current value of n
        directory = os.path.join(base_dir, f"bulk_cu_{n}x{n}x{n}")
        command = command_template.format(n)

        try:
            # Change the current working directory
            os.chdir(directory)
            print(f"Changed directory to {directory}")

            # Run the castep.serial command
            print(f"Running command: {command}")
            subprocess.run(command, shell=True, check=True)
        except Exception as e:
            print(f"An error occurred while processing {n}: {e}")


if __name__ == "__main__":
    run_castep_calculations(
        config.k_range
    )  # k_range should be defined in the config file
