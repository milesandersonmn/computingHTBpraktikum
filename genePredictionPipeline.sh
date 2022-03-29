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
geneModel=`basename -s .fa $TRANSCRIPT`
fileName=${genomeName}_${geneModel}

gth -genomic $GENOME -cdna $TRANSCRIPT -protein $PROTEIN \
-gff3out -o ${OUTPUT}/${fileName}.gff3 -startcodon -finalstopcodon -cdnaforward \
-skipalignmentout -v -exact -species maize -force

gffread -w ${OUTPUT}/${fileName}_exon.fa -g $GENOME ${OUTPUT}/${fileName}.gff3 &

gffread -C -x ${OUTPUT}/${fileName}_CDS.fa -g $GENOME ${OUTPUT}/${fileName}.gff3

seqkit translate ${OUTPUT}/${fileName}_CDS.fa > ${OUTPUT}/${fileName}_protein.fa

#Exon BLAST

blastn -subject ${OUTPUT}/${fileName}_exon.fa \
-query $TRANSCRIPT -outfmt "6 slen qlen sstart qstart send qend length pident\
  qcovus" > ${OUTPUT}/${fileName}_exon.blastn.txt &
  
blastn -subject $TRANSCRIPT \
-query ${OUTPUT}/${fileName}_exon.fa -outfmt "6 \
 qcovus" > ${OUTPUT}/${fileName}_exon.blastn.reverse.txt &
 
 #Protein BLAST
 
blastp -subject ${OUTPUT}/${fileName}_protein.fa \
-query $PROTEIN -outfmt "6 slen qlen sstart qstart send qend length pident \
qcovs" > ${OUTPUT}/${fileName}_protein.blastp.txt &
  
blastp -subject $PROTEIN \
-query ${OUTPUT}/${fileName}_protein.fa -outfmt "6 \
qcovs" > ${OUTPUT}/${fileName}_protein.blastp.reverse.txt &

#Concatenate BLAST and reverse BLAST

for i in `ls ${OUTPUT} | grep blastp.txt$`; do
sample=`basename -s .txt $i`
paste ${OUTPUT}/${sample}.txt ${OUTPUT}/${sample}.reverse.txt > \
${OUTPUT}/${sample}.combined.txt
done

for i in `ls ${OUTPUT} | grep blastn.txt$`; do
sample=`basename -s .txt $i`
paste ${OUTPUT}/${sample}.txt ${OUTPUT}/${sample}.reverse.txt > \
${OUTPUT}/${sample}.combined.txt
done


