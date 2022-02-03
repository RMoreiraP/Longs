#!/bin/bash

wd="$(pwd)"

invs="$(cut -f1 $wd/ListaRef.txt)"

# echo $invs

##  Bucle

# initBam="Fl_Mq_NMap_Merged.bam"
# initBam="SortedMapped_1.bam"
initBam="$wd/Fl_Mq30_NMap_HG002.bam"

invs="HsInv0370"

for inv in $invs
do
	$wd/changeTxt.sh ${inv}

	selInv=$(grep "$inv[[:space:]]" $wd/ListaRef.txt | cut -f2)
	grep "$inv[[:space:]]" $wd/ListaRef.txt | cut -f2 | sed 's/|/\t/g; s/://g; s/_/\t/g; s/-/\t/g; s/\./\t/g' | cut -f3,5,6 | sed 's/^[0]*//g' > $wd/$inv/seqInfo

	# echo $selInv

	## Extraigo ref

	samtools view -b $initBam $selInv > $wd/$inv/${inv}.bam
	samtools index $wd/$inv/$inv.bam > $wd/$inv/$inv.bam.bai

	# ##  Sacar reads completos
	# mkdir -p $wd/$inv/Resultados
	# > $wd/$inv/Resultados/ReadsA-D
	# 
	# initInv="$(cut -f1 $wd/$inv/CoordenadasRelativas/A | tr -d "[[:space:]]")"
	# finInv="$(cut -f2 $wd/$inv/CoordenadasRelativas/D | tr -d "[[:space:]]")"
	# 
	# samtools view $wd/$inv/$inv.bam --input-fmt-option "filter= pos < $initInv && endpos > $finInv" | cut -f1 > $wd/$inv/Resultados/ReadsA-D
	#
	# ## Eliminar reads que aparecen en ReadsA-D
	# 
	# samtools view -h $wd/$inv/$inv.bam | grep -vf $wd/$inv/Resultados/ReadsA-D | samtools view -bS -o $wd/$inv/remAD.bam
	# 
	# mv $wd/$inv/remAD.bam $wd/$inv/$inv.bam
	# rm $wd/$inv/$inv.bam.bai
	# samtools index $wd/$inv/$inv.bam > $wd/$inv/$inv.bam.bai

	## Mapear los reads en las sondas

	mkdir -p $wd/$inv/Resultados/SondasMapeadas
	letras="A B C D"

	> $wd/$inv/Resultados/MapInfo

	for sonda in $letras
	do
		initSonda="$(cut -f1 $wd/$inv/CoordenadasRelativas/$sonda | tr -d "[[:space:]]")"
		finSonda="$(cut -f2 $wd/$inv/CoordenadasRelativas/$sonda | tr -d "[[:space:]]")"

		samtools view $wd/$inv/$inv.bam -F 0x0010 --input-fmt-option "filter= pos <= $initSonda && endpos >= $finSonda" | cut -f1 > $wd/$inv/Resultados/SondasMapeadas/${sonda}_f
		samtools view $wd/$inv/$inv.bam -f 0x0010 --input-fmt-option "filter= pos <= $initSonda && endpos >= $finSonda" | cut -f1 > $wd/$inv/Resultados/SondasMapeadas/${sonda}_r

		> $wd/$inv/Resultados/SondasMapeadas/${sonda}

		for line in $(cat $wd/$inv/Resultados/SondasMapeadas/${sonda}_f)
		do
			echo -e "$line\t+" >> $wd/$inv/Resultados/SondasMapeadas/${sonda}
			echo -e "$line\t$sonda\t+" >> $wd/$inv/Resultados/MapInfo
		done

		for line in $(cat $wd/$inv/Resultados/SondasMapeadas/${sonda}_r)
		do
			echo -e "$line\t-" >> $wd/$inv/Resultados/SondasMapeadas/${sonda}
			echo -e "$line\t$sonda\t-" >> $wd/$inv/Resultados/MapInfo
		done

		rm $wd/$inv/Resultados/SondasMapeadas/${sonda}_*

	done



	## Tabla de T y F
	todosReads=$(cat $wd/$inv/Resultados/SondasMapeadas/* | cut -f1 | sort | uniq)
	> $wd/$inv/Resultados/TablaReadTF

	for mapRead in $todosReads
	do
		siA="F"
		siB="F"
		siC="F"
		siD="F"

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/A ;
		then
			siA="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/B ;
		then
			siB="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/C ;
		then
			siC="T"
		fi

		if grep -q $mapRead $wd/$inv/Resultados/SondasMapeadas/D ;
		then
			siD="T"
		fi

		echo -e "$mapRead\t$siA\t$siB\t$siC\t$siD" >> $wd/$inv/Resultados/TablaReadTF

	done

	#### En el siguiente comentario establezco que las "duplicadas" no lo son, que me coja todas las detectadas

	# for sonda in $letras; do cut -f1 $wd/$inv/Resultados/SondasMapeadas/$sonda | sort | uniq; done | sort | uniq -d > $wd/$inv/Resultados/readsDup
	for sonda in $letras; do cut -f1 $wd/$inv/Resultados/SondasMapeadas/$sonda | sort | uniq; done | sort | uniq > $wd/$inv/Resultados/readsDup

	## Filtrar para quedarnos solo con los reads importantes (quitar solo una sonda y reads solo BC)

	#### En el siguiente comentario establezco que no se filtran fuera los fragmentos BC

	# cat $wd/$inv/Resultados/TablaReadTF | grep -f $wd/$inv/Resultados/readsDup | awk '{if ($2 == "F" && $3 == "T" && $4 == "T" && $5 == "F"); else print}' > $wd/$inv/Resultados/TablaCandidatos
	cat $wd/$inv/Resultados/TablaReadTF | grep -f $wd/$inv/Resultados/readsDup > $wd/$inv/Resultados/TablaCandidatos
	
	cut -f1 $wd/$inv/Resultados/TablaCandidatos > $wd/$inv/Resultados/readsDup

	## Obtener BAM de reads que tienen varias sondas

	samtools view -H $wd/$inv/$inv.bam > $wd/$inv/header
	samtools view $wd/$inv/$inv.bam | grep -f $wd/$inv/Resultados/readsDup > $wd/$inv/reads

	cat $wd/$inv/header $wd/$inv/reads | samtools view -bS -o $wd/$inv/${inv}DupReads.bam
	rm $wd/$inv/header $wd/$inv/reads

	echo "Analyzing BLAST of $inv"
	$wd/anBlast.sh $inv


done

