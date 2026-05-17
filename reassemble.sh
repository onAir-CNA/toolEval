#!/bin/bash
IMAGE_NAME="jupyter-r-bioinfo:final"
PARTS_DIR="docker_image_parts"
OUTPUT_TAR="jupyter-r-bioinfo-final.tar.gz"

if docker images | grep -q "jupyter-r-bioinfo.*final"; then
    echo "Immagine $IMAGE_NAME già presente, skip riassemblaggio."
else
    echo "Riassemblaggio immagine da parti..."
    cat ${PARTS_DIR}/jupyter-r-bioinfo-final.tar.gz.part_* > ${OUTPUT_TAR}
    echo "Caricamento immagine Docker..."
    docker load < ${OUTPUT_TAR}
    rm -f ${OUTPUT_TAR}
    echo "Done! Immagine caricata come $IMAGE_NAME"
fi
