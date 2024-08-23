# pseudopotential_calculator

# Usage of the HPC

This package provides functions to help with
submitting jons to the HPC, however it is also sometimes to interact directly in the CLI using ssh

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
