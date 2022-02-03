#!/bin/bash

blastn -task blastn -query Ref.fasta -subject Ref.fasta -perc_identity 80 -outfmt "6 sstart send qstart qend qcovs" | awk '$1 != $3 && $2 != $4'
