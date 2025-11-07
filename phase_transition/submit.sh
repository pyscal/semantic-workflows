#!/bin/bash
#SBATCH --job-name=pdtest
#SBATCH --time=23:59:00
#SBATCH --partition=s.cmmg
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=3GB
#SBATCH --hint=nomultithread

python bcc.py
