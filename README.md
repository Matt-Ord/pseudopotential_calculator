# pseudopotential_calculator

## Interacting with the HPC

This package also contains functions to upload files to the HPC, and submit jobs remotely.

The behavior of these functions can be configured using environment variables

- `HPC_ACCOUNT`: The account used to submit jobs on the HPC
- `HPC_USERNAME` : The user used to sign in to the HPC
- `HPC_WORKSPACE` : The workspace used to upload files, the default being `/rds/user/{HPC_USERNAME}/hpc-work`

- `SKIP_DOWNLOAD` : The boolean variable indicating whether to fetch files from hpc, when the function in analyze.py is called

- `SKIP_UPLOAD` : The boolean variable indicating whether to fetch files to hpc, when the function in prepare.py is called

By default, you will be promted each time these variables are
By default, you will be promted each time these variables are
required, but by setting the environment variable this can be skipped.

If using VSCODE, place this in your .vscode/launch.json to automatically set the `HPC_ACCOUNT` and `HPC_USERNAME` variable

```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "Python Debugger: Current File",
      "type": "debugpy",
      "request": "launch",
      "program": "${file}",
      "env": {
        "HPC_USERNAME": "YOUR_USERNAME",
        "HPC_ACCOUNT": "YOUR_ACCOUNT"
      },
      "console": "integratedTerminal"
    }
  ]
}
```

It is also recommended that you setup ssh keys, in order to skip the password prompt.

```bash
ssh-keygen -t ed25519 -C "username@cam.ac.uk"
eval `ssh-agent -s`
ssh-add ~/.ssh/id_ed25519
chmod 600 ~/.ssh/id_ed25519
ssh-copy-id -i ~/.ssh/id_ed25519.pub username@login.hpc.cam.ac.uk
```

## Directly interacting with the HPC

It is also sometimes useful to interact directly with the HPC in the CLI using ssh

```shell
ssh username@login.hpc.cam.ac.uk
```

Once connected to the HPC, these HPC specific commands
may be useful

```shell
# To see all queue
squeue
# To see the queue for user `username`
squeue -u username
# To see the predicted start time
# NOTE: this does not take into account jobs ending early,
# or other users jumping the queue
squeue -u username --start
# To cancel a job
scancel job_id
# To view your balance
mybalance
```

The following bash commands are also generally useful

```shell
# Change directory
cd ./relative/directory
# list directory
ls
# list all files, along with permissions
ls -al
# Make a file excecutable
chmod u+x filename
# Remove all files in data directory
rm -r data
# make a directory at data/copper/bulk
mkdir -p data/copper/bulk
```

To upload files to the HPC, the scp command should be used inside the local terminal

```shell
# Note: -r copies all files recursively
scp -r ./local_dir/ username@login.hpc.cam.ac.uk:target_dir
```
