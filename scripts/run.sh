#!/usr/bin/bash

if [ $# -eq 0 ]; then
    echo "Uso: bash run.sh <modo>"
    echo "  bash run.sh ITS"
    echo "  bash run.sh concatenated"
    exit 1
fi

MODO=$1

if [ "$MODO" == "ITS" ]; then
    python scripts/pipeline.py nucleotide importante/ITS.txt
    python scripts/desenhar_arvore.py tree_ITS/raxml.raxml.support tree_ITS/aligned.mrbayes.nex.con.tre ITS_001

elif [ "$MODO" == "concatenated" ]; then
    mkdir -p mafft

    for ficheiro in importante/*.txt; do
        marcador=$(basename "$ficheiro" .txt)
        python scripts/fetch_sequences.py nucleotide "$ficheiro"
        python scripts/simple_fasta.py downloaded.fasta
        mafft simplified.fasta > mafft/${marcador}_mafft.fasta
    done

    python scripts/Seqconcact.py mafft/*_mafft.fasta
    python scripts/pipeline.py concatenated.fasta
    python scripts/desenhar_arvore_concatenated.py tree_concatenated/raxml.raxml.support tree_concatenated/concatenated.mrbayes.nex.con.tre concatenated_v001
else
    echo "Erro: modo '$MODO' não reconhecido. Usa 'ITS' ou 'concatenated'."
    exit 1
fi
