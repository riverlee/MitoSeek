#genome/
    # Mitochondrion fasta sequences
    #1) rCRS assembly
    #http://www.ncbi.nlm.nih.gov/nuccore/251831106

    #2) hg19 assembly
    #http://www.ncbi.nlm.nih.gov/nuccore/NC_001807.4?report=genbank

    #) refGene annotation file
    #description of refGene is here at https://cgwb.nci.nih.gov/cgi-bin/hgTables?hgsid=79902&hgta_doSchemaDb=hg18&hgta_doSchemaTable=refGene
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
    
    gunzip refGene.txt.gz
    
    #get chromosome 1-22, and X,Y
    
   perl -lane 'if($F[2]=~/chr\d+$|chrX|chrY/) {@start=split ",",$F[9]; @end=split ",",$F[10];  for($i=0;$i<$F[8];$i++) { print join "\t",$F[2],$start[$i],$end[$i],$F[12]."|".$F[1]}}' refGene.txt >refGeneExon_withChr.bed
  
   perl -lane 'if($F[2]=~/chr\d+$|chrX|chrY/) {$F[2]=~s/chr//g;@start=split ",",$F[9]; @end=split ",",$F[10];  for($i=0;$i<$F[8];$i++) { print join "\t",$F[2],$start[$i],$end[$i],$F[12]."|".$F[1]}}' refGene.txt >refGeneExon_withoutChr.bed
   
   #How to get the total bases from a bed file
  
  











