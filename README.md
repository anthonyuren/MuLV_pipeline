### MuLV insertion site mapping pipeline for Webster et al. 
https://doi.org/10.1101/157800 

These are the scripts used to map MuLV ligation mediated PCR reactions to the mouse genome and build insertion sites with estimations of relative clonality based on number of sheared fragments and number of reads.

Original scripts were written by Bruce Bolt and Barbara Iadarola. Scripts were revised, documented and beta tested by Marian Dore.

---
trim.sh

Adapter triming

---
pipeline-full.sh

Aligns reads to genome

---
countbylinecombined-withoriqc-R1.pl

Create a table with key information
Input: sorted,mapped and paired bam file
output: Mate1.txt (ReadID, chr, LTRpos, Orient, InsertSize, LTRPos, Read1, ReadLen, and if comb exists)  and Mate2.txt (chr, LTR, #of uniq LTRs, orient, tot # of LTRs) 

---
group.R

LTR-genome junctions that map with 10 bases and are oriented in the same direction are grouped and treated as one "best base"
Input: Mat1.txt
Output: TOTs.txt

---
inserts.R

Add best base and median base information and combine all into one table of all insertion sites
Input: TOTs.txt
Output:inserts.txt

---
QC scripts 0 to 5
Create a master table of all inserts, identify inserts that are identical between samples (and hence potentially PCR artifacts or cross contamination) and filter out these from the master table, in some cases keeping the most clonal/earliest cloned insert.




