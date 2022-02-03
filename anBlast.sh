#!/bin/bash

wd="$(pwd)"

inv="$1"

if [ $# -eq 0 ]
then
	echo "No arguments"
	exit
fi

mkdir -p $wd/$inv/Resultados/Distancias

samtools fastq $wd/$inv/${inv}DupReads.bam > $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq

> $wd/$inv/Resultados/sondaErronea

for dread in $(cat $wd/$inv/Resultados/readsDup)
do
	mkdir -p $wd/$inv/Resultados/Distancias/Reads/$dread

	numfila=$(cat -n $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq | grep $dread | cut -f1 | tr -d " ")
	tail -n+$numfila $wd/$inv/Resultados/Distancias/${inv}DupReads.fastq | head -4 > $wd/$inv/Resultados/Distancias/Reads/$dread/read.fastq

	## Bucle por cada sonda que existe

	for sonda in A B C D
	do
		mkdir $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord
		blastn -query $wd/$inv/Sondas/${sonda} -subject $wd/$inv/Resultados/Distancias/Reads/${dread}/read.fastq -task blastn -outfmt "6 sstart send slen sstrand" -qcov_hsp_perc 90 -sorthits 4 > $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda}

		if [ $(cat $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda} | wc -l) -eq 0 ]
		then
			rm $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda}
		fi

		if [ $(cat $wd/$inv/Resultados/Distancias/Reads/${dread}/BLASTCoord/${sonda} | wc -l) -gt 1 ]
		then
			echo -e "${dread}\t${sonda}" >> $wd/$inv/Resultados/sondaErronea
		fi

	done

done
