# pseudopotential_calculator

## Interacting with the HPC

This package also contains functions to upload files to the HPC, and submit jobs remotely.

By default, the username must be provided each time you log on, however
by setting the HPC_USERNAME environment variable this can be skipped.

If using VSCODE, place this in your .vscode/launch.json to automatically set the username variable

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
        "HPC_USERNAME": "YOUR_USERNAME"
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
