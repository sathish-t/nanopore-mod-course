# Go-related path stuff, delete if not relevant for you.
export GOPATH=${HOME}/go
export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin

# using singularity images for packages we use for the course
export sgl_base_path=${HOME}/nanomod_course_singularity_img_collection

bedtools() {  singularity exec "$sgl_base_path/software.img" bedtools "$@"; }
export -f bedtools;

DNAscent() { singularity exec "$sgl_base_path/dnascent.sif" /app/DNAscent/bin/DNAscent "$@"; }
export -f DNAscent;

dorado() { singularity exec "$sgl_base_path/software.img" dorado "$@"; }
export -f dorado;

minimap2() { singularity exec "$sgl_base_path/software.img" minimap2 "$@"; }
export -f minimap2;

modbamtools() { singularity exec "$sgl_base_path/software.img" modbamtools "$@"; }
export -f modbamtools;

modkit() { singularity exec "$sgl_base_path/software.img" modkit "$@"; }
export -f modkit;

nanalogue() { singularity exec "$sgl_base_path/software.img" nanalogue "$@"; }
export -f nanalogue;

pycoQC() { singularity exec "$sgl_base_path/software.img" pycoQC "$@"; }
export -f pycoQC;

samtools() { singularity exec "$sgl_base_path/software.img" samtools "$@"; }
export -f samtools;

aws() { singularity exec "$sgl_base_path/software.img" aws "$@"; }
export -f aws;

pod5() { singularity exec "$sgl_base_path/software.img" pod5 "$@"; }
export -f pod5;

nanalogue-gui() { cd ~/nanomod_course_software/nanalogue-gui && npm start; }
export -f nanalogue-gui;

# Not sure why following line is needed
export LD_PRELOAD=

# Some NVM related stuff for nanalogue-gui
export NVM_DIR="$HOME/.nvm"
[ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"  # This loads nvm
[ -s "$NVM_DIR/bash_completion" ] && \. "$NVM_DIR/bash_completion"  # This loads nvm bash_completion
