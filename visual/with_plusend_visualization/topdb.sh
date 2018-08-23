#!/bin/tcsh
#SBATCH --mail-user=qni@umd.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 5:00:00
#SBATCH --mem 5120
rm -f topdb.out
module load matlab
matlab -nodisplay -nosplash < runfile.m > ./topdb.out

