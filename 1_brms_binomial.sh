#!/bin/bash

#Path

cd /home/nbogdanovic/WNV_prevalence

# Loading program
module load R-bundle-CRAN
# module load git/2.45.1

# Pull commits
# git pull

# Sending the job
R CMD BATCH --vanilla 1_brms_binomial.R 1_brms_binomial_random.Rout

# Commit and push the files
# git add --all
# git commit -m "job from cluster: Preparing_landcover_grid_shp"
# git pull
# git push


# run using:
# sbatch -p ceab --ntasks=4 --mem=32G --mail-type=BEGIN,END,FAIL --mail-#user=nina.bogdanovic@ceab.csic.es 1_brms_binomial.sh


