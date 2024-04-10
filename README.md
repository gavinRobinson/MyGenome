# MyGenome
Analyses for ABT480/CS485G genome assembly

## 1. Analysis of sequence quality
The F1 and R1 sequence data sets were analyzed using FASTQC:
```bash
ssh -Y ghro223@ghro223.cs.uky.edu
cd MyGenome
fastqc&
```
Load F1 and R1 data into GUI interface.
Take screen shots of output files: 

![Screenshot(111).png](data/U247_1_paired.png)
![Screenshot(111).png](data/U247_2_paired.png)

## 2. Assembly using MCC
Logging into MCC and scp sequence data
```bash
ssh ghro223@mc.uky.edu
cd..
cd..
cd project/farman_s24cs485g/
mkdir ghro223
cd ghro223
scp ghro223@ghro223.cs.uky.edu:/home/ghro223/MyGenome/U247_1_paired.fastq .
scp ghro223@ghro223.cs.uky.edu:/home/ghro223/MyGenome/U247_1_unpaired.fastq .
scp ghro223@ghro223.cs.uky.edu:/home/ghro223/MyGenome/U247_2_paired.fastq .
scp ghro223@ghro223.cs.uky.edu:/home/ghro223/MyGenome/U247_2_unpaired.fastq .
cp ../SLURM_SCRIPTS/velvetoptimiser_noclean.sh .
vim velvetoptimiser_noclean.sh
```
I then altered the velvetoptimiser_noclean.sh script to include my email and some new header data, and then went on to run the velvetoptimiser_noclean.sh with a step size of 10
```bash
sbatch velvetoptimiser_noclean.sh U247 61 131 10
```
After the run was complete, I analyzed the new assembly data via command
```bash
tail -50 U247/velvet_U247_61_131_10_noclean/19-03-2024-15-01-19_Logfile.txt
``` 
![Screenshot(121).png](data/Screenshot(121).png) 
I then went on to run the velvetoptimiser_noclean.sh script with a step size of 2
```bash
sbatch velvetoptimiser_noclean.sh U247 110 130 2
```
