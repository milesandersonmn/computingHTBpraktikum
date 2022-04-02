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

gffread -w ${OUTPUT}/${fileName}_exon.fa -g $GENOME ${OUTPUT}/${fileName}.gff3

gffread -C -x ${OUTPUT}/${fileName}_CDS.fa -g $GENOME ${OUTPUT}/${fileName}.gff3

seqkit translate ${OUTPUT}/${fileName}_CDS.fa > ${OUTPUT}/${fileName}_protein.fa

#Exon BLAST

blastn -subject ${OUTPUT}/${fileName}_exon.fa \
-query $TRANSCRIPT -outfmt "6 sseqid slen qlen sstart qstart send qend length pident\
  qcovus" > ${OUTPUT}/${fileName}_exon.blastn.txt &&
  
blastn -subject $TRANSCRIPT \
-query ${OUTPUT}/${fileName}_exon.fa -outfmt "6 \
 qcovus" > ${OUTPUT}/${fileName}_exon.blastn.reverse.txt &&
 
 #Protein BLAST
 
blastp -subject ${OUTPUT}/${fileName}_protein.fa \
-query $PROTEIN -outfmt "6 sseqid slen qlen sstart qstart send qend length pident \
qcovs" > ${OUTPUT}/${fileName}_protein.blastp.txt &&
  
blastp -subject $PROTEIN \
-query ${OUTPUT}/${fileName}_protein.fa -outfmt "6 \
qcovs" > ${OUTPUT}/${fileName}_protein.blastp.reverse.txt &&

#Concatenate BLAST and reverse BLAST

paste ${OUTPUT}/${fileName}_protein.blastp.txt ${OUTPUT}/${fileName}_protein.blastp.reverse.txt > \
${OUTPUT}/${fileName}_protein.blastp.combined.txt

paste ${OUTPUT}/${fileName}_exon.blastn.txt ${OUTPUT}/${fileName}_exon.blastn.reverse.txt > \
${OUTPUT}/${fileName}_exon.blastn.combined.txt

#Count number of mRNA transcripts and proteins output

transcriptCount=`cat ${OUTPUT}/${fileName}_exon.blastn.combined.txt | awk '{print $1}' | uniq | wc -l`

proteinCount=`cat ${OUTPUT}/${fileName}_protein.blastp.combined.txt | awk '{print $1}' | uniq | wc -l`

echo "Number of predicted mRNA transcripts: ${transcriptCount} \n Number of predicted proteins: ${proteinCount}"



