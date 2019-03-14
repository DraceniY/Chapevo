
#!/bin/bash

	for FastaFile in *.fna ; do
		echo 'WORKING ON FILE:' $FastaFile
		echo 'BlastN'

		newFASTA=$(echo $FastaFile | cut -f 1 -d '.')
		blastn -query $FastaFile -db "/media/yasmine/DATA2/RnaR_18S/Blast/18S_HomoS" -out $newFasta"_18S.out" -evalue 0.000001 -outfmt 6

	done
	mkdir output
	mv *.out output

#Clean up







