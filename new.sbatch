#!/bin/bash 
#dhrf.sbatch, a sample slurm job  
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=sim 
#################  
#a file for job output, you can check job progress
#SBATCH --output=sim.out
#################
# a file for errors from the job
#SBATCH --error=sim.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH --time=2880:00
#################
#quality of service; think of it as job priority
#################
#number of nodes you are requesting
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=enf@stanford.edu
#################
##SBATCH -n 2 -p rondror --gres=gpu:2 --qos rondror
#SBATCH -n 2 -p gpu --qos=gpu --gres=gpu:2 --constraint="titanblack"
#################
#now run normal batch commands

#SBATCH --dependency=singleton

if [[ ! -f Min_3.rst ]]; then
    ml load amber/14-intel
    srun pmemd -O -i Min_1.in -o Min_1.out -p system.prmtop -c system.inpcrd -r Min_1.rst -ref system.inpcrd -x Min_1.nc
    srun pmemd -O -i Min_2.in -o Min_2.out -p system.prmtop -c Min_1.rst -r Min_2.rst -ref Min_1.rst -x Min_2.nc
    srun pmemd -O -i Min_3.in -o Min_3.out -p system.prmtop -c Min_2.rst -r Min_3.rst -ref Min_2.rst -x Min_3.nc
    ml load amber/14-cuda
    srun pmemd.cuda.MPI -O -i Heat_1.in -o Heat_1.out -p system.prmtop -c Min_3.rst -r Heat_1.rst -ref Min_3.rst -x Heat_1.nc
    srun pmemd.cuda.MPI -O -i Heat_2.in -o Heat_2.out -p system.prmtop -c Heat_1.rst -r Heat_2.rst -ref Min_3.rst -x Heat_2.nc
    srun pmemd.cuda.MPI -O -i Eq_1.in -o Eq_1.out -p system.prmtop -c Heat_2.rst -r Eq_1.rst -ref Min_3.rst -x Eq_1.nc
    srun pmemd.cuda.MPI -O -i Eq_2.in -o Eq_2.out -p system.prmtop -c Eq_1.rst -r Eq_2.rst -ref Min_3.rst -x Eq_2.nc
    srun pmemd.cuda.MPI -O -i Eq_3.in -o Eq_3.out -p system.prmtop -c Eq_2.rst -r Eq_3.rst -ref Min_3.rst -x Eq_3.nc
    srun pmemd.cuda.MPI -O -i Eq_4.in -o Eq_4.out -p system.prmtop -c Eq_3.rst -r Eq_4.rst -ref Min_3.rst -x Eq_4.nc
    srun pmemd.cuda.MPI -O -i Eq_5.in -o Eq_5.out -p system.prmtop -c Eq_4.rst -r Eq_5.rst -ref Min_3.rst -x Eq_5.nc
    srun pmemd.cuda.MPI -O -i Eq_6.in -o Eq_6.out -p system.prmtop -c Eq_5.rst -r Eq_6.rst -ref Min_3.rst -x Eq_6.nc
    srun pmemd.cuda.MPI -O -i Prod.in -o Prod_1.out -p system.prmtop -c Eq_6.rst -r Prod_1.rst -ref Min_3.rst -x Prod_1.nc
elif [[ ! -f Eq_6.rst ]]; then  
    ml load amber/14-cuda
    srun pmemd.cuda.MPI -O -i Heat_1.in -o Heat_1.out -p system.prmtop -c Min_3.rst -r Heat_1.rst -ref Min_3.rst -x Heat_1.nc
    srun pmemd.cuda.MPI -O -i Heat_2.in -o Heat_2.out -p system.prmtop -c Heat_1.rst -r Heat_2.rst -ref Min_3.rst -x Heat_2.nc
    srun pmemd.cuda.MPI -O -i Eq_1.in -o Eq_1.out -p system.prmtop -c Heat_2.rst -r Eq_1.rst -ref Min_3.rst -x Eq_1.nc
    srun pmemd.cuda.MPI -O -i Eq_2.in -o Eq_2.out -p system.prmtop -c Eq_1.rst -r Eq_2.rst -ref Min_3.rst -x Eq_2.nc
    srun pmemd.cuda.MPI -O -i Eq_3.in -o Eq_3.out -p system.prmtop -c Eq_2.rst -r Eq_3.rst -ref Min_3.rst -x Eq_3.nc
    srun pmemd.cuda.MPI -O -i Eq_4.in -o Eq_4.out -p system.prmtop -c Eq_3.rst -r Eq_4.rst -ref Min_3.rst -x Eq_4.nc
    srun pmemd.cuda.MPI -O -i Eq_5.in -o Eq_5.out -p system.prmtop -c Eq_4.rst -r Eq_5.rst -ref Min_3.rst -x Eq_5.nc 
    srun pmemd.cuda.MPI -O -i Eq_6.in -o Eq_6.out -p system.prmtop -c Eq_5.rst -r Eq_6.rst -ref Min_3.rst -x Eq_6.nc
    srun pmemd.cuda.MPI -O -i Prod.in -o Prod_1.out -p system.prmtop -c Eq_6.rst -r Prod_1.rst -ref Min_3.rst -x Prod_1.nc
elif [[ -f Prod_1.rst ]] || [[ -f Prod.rst ]]; then 
    last=$(ls -1 Prod_*.out | sed 's/.*_\([0-9]\+\).*/\1/' | sort -n | tail -n1)
    ml load amber/14-cuda
    srun pmemd.cuda.MPI -O -i Prod.in -o Prod_$((last+1)).out -p system.prmtop -c Prod_$last.rst -r Prod_$((last+1)).rst -ref Prod_$last.rst -x Prod_$((last+1)).nc
fi
    echo "else condition found"
