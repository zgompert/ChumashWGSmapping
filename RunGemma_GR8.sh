#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=220G
#SBATCH --account=gompert
#SBATCH --partition=gompert-grn
#SBATCH --qos=gompert-grn
#SBATCH --job-name=gemma
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load gemma
# 0.95a
cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gemma_rev

gfile="hsub_mlpntest_GR8_tchum_gr8_green.geno"
pfile="ph_tchumGR8.txt"
oseed="o_GR8_"

#gfile="hsub_mlpntest_tchum_red.geno"
#pfile="ph_tchum_red.txt"
#oseed="o_red_"
#gfile="hsub_mlpntest_tchum_green.geno"
#pfile="ph_tchum_green.txt"
#oseed="o_green_"

## RG 
of="{$oseed}RG"
gemma -g  $gfile -p $pfile -gk 1 -o $of -maf 0
gemma -g  $gfile -p $pfile -k output/$of.cXX.txt -lmm 4 -n 1 -o $of -maf 0

## GB 
of="{$oseed}GB"
gemma -g  $gfile -p $pfile -gk 1 -o $of -maf 0
gemma -g  $gfile -p $pfile -k output/$of.cXX.txt -lmm 4 -n 2 -o $of -maf 0

