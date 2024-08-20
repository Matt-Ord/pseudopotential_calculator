import os

from ase.build import bulk
from ase.calculators.castep import Castep

# Define bulk copper structure
bulk_copper = bulk("Cu", "fcc", a=3.8)  # FCC copper, a = 3.615
for k in range(1, 11):
    directory = f"data/bulk_cu_{k}x{k}x{k}"
    label = f"bulk_cu_{k}x{k}x{k}"
    slurm_script_path = os.path.join(directory, f"submit_job.sh")

    calculation = Castep(
        directory=directory,
        label=label,
        castep_command="castep.serial",
        check_castep_version=True,
        keyword_tolerance=0,
    )
    calculation.xc_functional = "PBE"
    calculation.cut_off_energy = 340  # in eV
    calculation.cell.kpoint_mp_grid = f"{k} {k} {k}"
    # Controlling the angles within the cell so that it's FCC all the time during simulation
    calculation.cell.cell_constraints = "1 1 1\n0 0 0"
    calculation.task = "GeometryOptimization"

    calculation._track_output = False # type: ignore : bad library type
    calculation._try_reuse = True
    calculation._pedantic = True
    calculation._rename_existing_dir = False
    calculation.param.reuse = True
    calculation.param.num_dump_cycles = 0
    calculation._export_settings = True
    # Attach the calculator to the atoms object
    calculation.set_atoms(bulk_copper)
    bulk_copper.set_calculator(calculation)

    if calculation.dryrun_ok():
        print("ok")

        calculation.prepare_input_files()
        # Write the SLURM submit script
        with open(slurm_script_path, "w") as f:
            f.write(
                """
    #! sbatch directives begin here ###############################
    #! Name of the job:
    #SBATCH -J copper_bulk
    #! Which project should be charged:
    #SBATCH -A CHANGEME
    #SBATCH -p cclake
    #! How many whole nodes should be allocated?
    #SBATCH --nodes=1
    #! How many (MPI) tasks will there be in total? (<= nodes*56)
    #! The Cascade Lake (cclake) nodes have 56 CPUs (cores) each and
    #! 3410 MiB of memory per CPU.
    #SBATCH --ntasks=56
    #! How much wallclock time will be required?
    #SBATCH --time=02:00:00
    #! What types of email messages do you wish to receive?
    #SBATCH --mail-type=NONE
    #! Uncomment this to prevent the job from being requeued (e.g. if
    #! interrupted by node failure or system downtime):
    ##SBATCH --no-requeue

    #! sbatch directives end here (put any additional directives above this line)

    #! Notes:
    #! Charging is determined by cpu number*walltime.
    #! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
    #! usually equates to the number of MPI tasks launched. Reduce this from nodes*56 if
    #! demanded by memory requirements, or if OMP_NUM_THREADS>1.
    #! Each task is allocated 1 CPU by default, and each CPU is allocated 3420 MiB
    #! of memory. If this is insufficient, also specify
    #! --cpus-per-task and/or --mem (the latter specifies MiB per node).

    #! Number of nodes and tasks per node allocated by SLURM (do not change):
    numnodes=$SLURM_JOB_NUM_NODES
    numtasks=$SLURM_NTASKS
    mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\\([0-9][0-9]*\\).*$/\1/')
    #! ############################################################
    #! Modify the settings below to specify the application's environment, location
    #! and launch method:

    #! Optionally modify the environment seen by the application
    #! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
    . /etc/profile.d/modules.sh                # Leave this line (enables the module command)
    module purge                               # Removes all modules still loaded
    module load rhel8/default-ccl              # REQUIRED - loads the basic environment

    #! Insert additional module load commands after this line if needed:
    """,
            )
# write a submit_all_script

with open(data, "w") as f:
            f.write(
    f"""
    for i in range(1,11):
        cd ./bulk_cu_{k}x{k}x{k}
        sbatch- submit_job.sh
        cd ..

    """,
            )
    else:
        msg = "dryrun failed"
        raise Exception(msg)
