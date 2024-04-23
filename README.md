# MyGenome (U247)
Analyses for ABT480/CS485G genome assembly

## 1. Analysis of U247 sequence quality
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

## 2. Assembly of U247 using MCC
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
After the run was complete, I analyzed the new assembly data via command
```bash
tail -50 U247/velvet_U247_110_130_2_noclean/21-03-2024-14-41-21_Logfile.txt
```
![Screenshot(127).png](data/Screenshot(127).png)
## 3. Blasting U247
Inside of my blast directory on my VM I ran the following
```bash
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn0 -evalue 1e-20 -outfmt 0
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn1 -evalue 1e-20 -outfmt 1
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn2 -evalue 1e-20 -outfmt 2
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn3 -evalue 1e-20 -outfmt 3
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn4 -evalue 1e-20 -outfmt 4
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn5 -evalue 1e-20 -outfmt 5
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn6 -evalue 1e-20 -outfmt 6
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn7 -evalue 1e-20 -outfmt 7
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn8 -evalue 1e-20 -outfmt 8
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn9 -evalue 1e-20 -outfmt 9
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn10 -evalue 1e-20 -outfmt 10
blastn -subject U247.fasta-query MoRepeats.fasta -out MoRepeats.U247.BLASTn11 -evalue 1e-20 -outfmt 11
```
I then used scp to transfer SequenceLengths.pl from the farman lab computer desktop to my VM and ran the following in my blast directory
```bash
perl SequenceLengths.pl U247.fasta | sort -k2n
```
![Screenshot(130).png](data/Screenshot(130).png)
From this I learned that my longest contig was U247_contig2655.
I then ran the following on my n6 output in an attempt to see if my genome contained any MAGGY gene
```bash
grep MAGGY MoRepeats.U247.BLASTn6
```
And did not get any matches, I then ran the same but with the Pot2 gene
```bash
grep Pot2 MoRepeats.U247.BLASTn6
grep Pot2 MoRepeats.U247.BLASTn6 | wc -l
```
And got 85 matches
![Screenshot(129).png](data/Screenshot(129).png)
I then ran the following because of the assignment and yielded no output
```blast
grep Pot2 MoRepeats.U247.BLASTn6 | awk '$4 >= 5638*0.9'
grep Pot2 MoRepeats.U247.BLASTn6 | awk '$2 ~ /contig2655/' | awk '$9 > 2000000 && $9 < 3000000' | sort -k9n
```
I then ran the following to remove contigs from my assembly that are less than 200 base pairs in length
```bash
perl CullShortContigs.pl U247_nh.fasta
perl SequenceLengths.pl U247_final.fasta
```
The following code was to determine which contigs in the resulting assembly correspond to the mitochondrial genome (NCBI needs us to identify those contigs).
BLAST the MoMitochondrion.fasta sequence against your final genome assembly using output format 6 with specific column selections
```bash
blastn -query MoMitochondrion.fasta -subject U247_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' -out MoMitochondrion.U247.BLAST
```
Then exporting a list of contigs that mostly comprise mitochondrial sequences (this will be uploaded to NCBI along with our genome assembly):
```bash
awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' MoMitochondrion.U247.BLAST > U247_mitochondrion.csv
```
This command yielded an output [.csv file called U247_mitochondrion.csv](data/U247_mitochondrion.csv) which is to be submitted, along with my U247_final.fasta, for the NCBI submission. I then copied the B71v2sh_masked.fasta genome from (Farman Lab machine Desktop) to my blast directory inside my VM and ran the following 
```bash
blastn -query B71v2sh_masked.fasta -subject U247_Final.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.U247.BLAST
```
I then created a directory named MyGenome_BLAST inside of my MCC ghro223 directory and moved my BLAST output file into that directory. I also copied the B71v2sh.U247.BLAST file into the CLASS_BLASTS directory (/project/farman_s24cs485g/BLAST). I then copied the CallVariants.sh script from farman_s24cs485g/SCRIPTs into the ghro223 directory in farman_s24cs485g and ran the following SLURM script on the MCC
```bash
sbatch CallVariants.sh U247_BLAST
```
## 4. Busco U247
## 5. Gene Prediction
