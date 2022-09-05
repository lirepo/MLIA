#!/bin/bash





#run path
RUN_PATH=/share/home/pi_zgli/zgli/mito_eq/1.assembly

# pacbio clr read
PACBIO_FQ=/share/home/pi_zgli/zgli/mito_eq/data/Eq.pb.fastq.gz

# mito ref genome
REF_MITO_FA=Erysiphe_pisi.fasta

#cpu number
CPU=20

#filter cutoff
HIT_MATCH_LEN_CUTOFF=1000
MIN_READ_LEN_CUTOFF=2000


# mito genome size
MITO_SIZE=200000


#iter number
ITER_NUM=3

# init read id and sequence for the 1st loop
firstReadId=4.readHitId.iter0.union.list


######################################################################################################

cd ${RUN_PATH}
touch ${firstReadId}

for i in $(seq 1 ${ITER_NUM})
do


# 1.run minimap2 searching against ref mito genome 
COMMAND="minimap2  -t ${CPU}  -x map-pb ${REF_MITO_FA} ${PACBIO_FQ} -o 1.minimap.iter${i}.out.paf"
echo -ne "CUSTOME_RUN:	iter${i} 1\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}


# 2. filter aligned hit
COMMAND="awk '\$10> ${HIT_MATCH_LEN_CUTOFF} && \$2 > ${MIN_READ_LEN_CUTOFF}'  1.minimap.iter${i}.out.paf  >  2.minimap_filter.iter${i}.paf"
echo -ne "CUSTOME_RUN:	iter${i} 2\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}


# 3. get read ids that hit the genome
COMMAND="cut -f 1 2.minimap_filter.iter${i}.paf | sort -u > 3.readHitId.iter${i}.list"
echo -ne "CUSTOME_RUN:	iter${i} 3\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}


# 4. combine readId with readId previous identified
previous_loop_i=$((${i}-1))
COMMAND="cat 3.readHitId.iter${i}.list  4.readHitId.iter${previous_loop_i}.union.list >4.readHitId.iter${i}.union.list"
echo -ne "CUSTOME_RUN:	iter${i} 4\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}


# uniq read ids
COMMAND="sort -u 4.readHitId.iter${i}.union.list > 4.readHitId.iter${i}.unionUniq.list"
echo -ne "CUSTOME_RUN:	iter${i} 4\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}




# 5. get read sequence from Pacbio FQ file
COMMAND="seqtk subseq ${PACBIO_FQ} 4.readHitId.iter${i}.unionUniq.list > 5.readHitSeq.iter${i}.union.fq"
echo -ne "CUSTOME_RUN:	iter${i} 5\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}




# 6. run canu to assembly 
COMMAND="canu -p mito -d 6.canu.iter${i} genomeSize=${MITO_SIZE} -pacbio 5.readHitSeq.iter${i}.union.fq"
echo -ne "CUSTOME_RUN:	iter${i} 6\t"
echo -ne `date +%Y-%m-%d,%H:%m:%s`
echo -ne "\t"
echo ${COMMAND}
#eval ${COMMAND}


# 7. set the new ref mito genome  
REF_MITO_FA=6.canu.iter${i}/mito.contigs.fasta


done
