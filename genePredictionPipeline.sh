declare option
declare OPTARG
OPTIND=1
GENOME=''
TRANSCRIPT=''
PROTEIN=''
OUTPUT=''

print_help() { echo "Usage: sh Annotation.sh -g path/to/genome/fasta -t path/to/transcript/file -p path/to/protein/file \
-o path/to/output/directory" 1>&2; exit 1; }

while getopts g:t:p:o:h option
do
case "${option}"
in
g) GENOME=${OPTARG};;
t) TRANSCRIPT=${OPTARG};;
p) PROTEIN=${OPTARG};;
o) OUTPUT=${OPTARG};;
h) print_help; exit 2;;
*) print_help; exit 2;;
esac
done

genomeName=`basename -s .fa $GENOME`
geneModel=`basename -s ,fa $TRANSCRIPT`
fileName=${genomeName}_${geneModel}

gth -genomic $GENOME -cdna $TRANSCRIPT -protein $PROTEIN \
-gff3out -o $OUTPUT/$fileName.gff3 -startcodon -finalstopcodon -cdnaforward \
-skipalignmentout -v -exact -species maize -force
