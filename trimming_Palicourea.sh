for infile in /media/ambedoya/extradrive1/Ana/Palicourea/raw_data/*_1.fq.gz
	do
	base=$(basename ${infile} _1.fq.gz)
	trimmomatic PE -threads 10 ${infile} ${base}_2.fq.gz ${base}_1.trim.fq.gz ${base}_1un.trim.fq.gz ${base}_2.trim.fq.gz ${base}_2un.trim.gz ILLUMINACLIP:novogene_adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	done

