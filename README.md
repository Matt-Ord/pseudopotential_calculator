# pseudopotential_calculator

## Interacting with the HPC

This package also contains functions to upload files to the HPC, and submit jobs remotely.

By default, the username must be provided each time you log on, however
by setting the HPC_USERNAME environment variable this can be skipped.

If using VSCODE, place this in your .vscode/launch.json to automatically set the HPC_ACCOUNT variable

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
        "HPC_ACCOUNT": "YOUR_USERNAME"
      },
      "console": "integratedTerminal"
    }
  ]
}
```

The package also supports custom HPC workspace directories, which can be modified by setting
the HPC_WORKSPACE variable. This defaults to `/rds/user/{HPC_USERNAME}/hpc-work`.

It is also recommended that you setup ssh keys, in order to skip the password prompt.

```bash
ssh-keygen -t ed25519 -C "username@cam.ac.uk"
eval `ssh-agent -s`
ssh-add ~/.ssh/id_ed25519
chmod 600 ~/.ssh/id_ed25519
ssh-copy-id -i ~/.ssh/id_ed25519.pub username@login.hpc.cam.ac.uk
```

# Directly interacting with the HPC

It is also sometimes useful to interact directly with the HPC in the CLI using ssh

```shell
ssh username@login.hpc.cam.ac.uk
```

Once connected to the HPC, these HPC specific commands
may be useful

```shell
# To see all queue
squeue
# To see the queue for user username
squeue -u username
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
```

To upload files to the HPC, the scp command should be used inide the local terminal

```shell
# Note: -r copies all files recursively
scp -r ./local_dir username@login.hpc.cam.ac.uk:target_dir
```
