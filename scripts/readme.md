```shell
# To login
$ ssh ll720@login.hpc.cam.ac.uk
# To get current path
$ echo "$(pwd -P)"
# The module setting inside submit script
module load rhel7/default-peta4 castep          # REQUIRED - loads the basic environment
# To upload repository. Do it in local terminal. Not in the remote one.
$ scp -r ~/local_dir user@host.com:/var/www/html/target_dir
$ scp -r ./data ll720@login.hpc.cam.ac.uk:/rds/user/ll720/hpc-work
# To remove repository
$ rm -r data
# To grant execution permission
$ chmod +x file_name.sh
# To run a bash script
$ ./submit_all.sh
# To see all queue
$ squeue
# To see my queue
$ squeue -u ll720
# To remove a job
$ scancel job_id

# To check my balance
$ my  balance
# To retreive file from remote
$ scp -r ll720@login.hpc.cam.ac.uk:rds/hpc-work/data/bulk_cu_1x1x1 ./data

# For Permission
    # To review permission
    $ ls -alh
    # To edit permission
    $ chmod u+x filename
    # To edit permission for all files within a directory
    $ chmod -R +x filename
/workspaces/pseudopotential_calculator/data
```
