#!/bin/bash
for i in *R1.fastq.gz
 do
 base=${i%_*.gz}
 
 #align all reads with MompS1 including flanking region
 #bwa mem Ref_Paris_mompS_1_flank2.fasta ${base}_R1.fastq.gz ${base}_R2.fastq.gz > ${base}-pe.sam
 bwa mem -t 16 Ref_Paris_mompS_1_flank2.fasta ${base}_R1.fastq.gz ${base}_R2.fastq.gz | samtools sort -@16 -O BAM -o ${base}_aln-pe.sorted.bam 
 
 #sort and index sam file
 #samtools view -bS ${base}-pe.sam > ${base}_aln-pe.bam
 #amtools sort ${base}_aln-pe.bam -o ${base}_aln-pe.sorted.bam
 #samtools index ${base}_aln-pe.sorted.bam
 parallel samtools index {} ::: *.sorted.bam
 
 #select only reads which start or end in anchoring region by ID
 samtools view -f 0x2 -L mompS1-positions.bed ${base}_aln-pe.sorted.bam | cut -f1 > ${base}_IDs.txt
 java -jar picard.jar FilterSamReads I=${base}_aln-pe.sorted.bam O=${base}-pe_final.bam READ_LIST_FILE=${base}_IDs.txt FILTER=includeReadList
 samtools sort ${base}-pe_final.bam -o ${base}-pe_final.sorted.bam
 #samtools index ${base}-pe_final.sorted.bam
 parallel samtools index {} ::: *.sorted.bam
 
 #Obtain consensus sequence from retained aligned reads
 samtools mpileup -uf Ref_Paris_mompS_1_flank2.fasta ${base}-pe_final.sorted.bam | bcftools call -c --ploidy 1 -Ov -o ${base}-calls.vcf
 bgzip ${base}-calls.vcf
 tabix -f *calls.vcf.gz
 bcftools consensus -f Ref_Paris_mompS_1_flank2.fasta ${base}-calls.vcf.gz > ${base}-cns.fa
 
 #Clean up
 rm ${base}-pe.sam
 rm ${base}_aln-pe.sorted.bam*
 rm ${base}_aln-pe.bam
 rm ${base}*reads*
done

#Use mlst to type mompS1

