# pacbio clr read
PACBIO_FQ = "pb.fq.gz"

# mito ref genome
REF_MITO_FA ="Erysiphe_pisi.fasta"

#cpu number
CPU = 20

# mito genome size
MITO_SIZE = 200000

#iter number
ITER_NUM = 5


rule all:
    input:
        "iter_"+ str(ITER_NUM) +".step_6.Canu/mito.contigs.fasta"


def looper_input(wildcards):
    n = int(wildcards.n)
    if n == 1:
        return REF_MITO_FA
    elif n > 1:
        return "iter_%d.step_6.Canu/mito.contigs.fasta" % (n - 1)
    else:
        raise ValueError("Loop numbers must be 1 or greater: received %d" % n)

wildcard_constraints:
    n="[0-9]+"

# 1.run minimap2 searching against ref mito genome
rule minimap2_search:
    input:
        looper_input
    output:
        "iter_{n}.step_1.minimap2/map.paf",
    shell:
        "minimap2  -t {CPU}  -x map-pb {input} {PACBIO_FQ} -o {output}"
        "mkdir -p iter_0.step_3.hitId"
        "touch iter_0.step_3.hitId/readHitId.list"

# 2. filter aligned hit
rule filterHit:
    input:
        "iter_{n}.step_1.minimap2/map.paf"
    output:
        "iter_{n}.step_2.filterHit/map_filter.paf"
    params:
        HIT_MATCH_LEN_CUTOFF = 1000,
        MIN_READ_LEN_CUTOFF = 2000
    shell:
        "awk '\$10> {params.HIT_MATCH_LEN_CUTOFF} && \$2 > {params.MIN_READ_LEN_CUTOFF}'  {input}  >  {output}"


# 3. get read ids that hit the genome
rule getReadId:
    input:
        "iter_{n}.step_2.filterHit/map_filter.paf"
    output:
        "iter_{n}.step_3.hitId/readHitId.list"
    shell:
        "cut -f 1 {input} | sort -u > {output}"

# 4. combine readId with readId previous identified
def step4_previousIter(wildcards):
    n = int(wildcards.n)
    return "iter_%d.step_3.hitId/readHitId.list" % (n-1)

rule readIdUnion:
    input:
        "iter_{n}.step_3.hitId/readHitId.list",
    output:
        "iter_{n}.step_4.unionReadId/unionUniq.list"
    params:
        readId_previsou_iter = step4_previousIter
    shell:
        "cat {input} {params.readId_previsou_iter} | sort -u > {output}"


# 5. get read sequence from Pacbio FQ file
rule getReadSeq:
    input:
        "iter_{n}.step_4.unionReadId/unionUniq.list"
    output:
        "iter_{n}.step_5.readHitSeq/union.fq"
    shell:
        "seqtk subseq {PACBIO_FQ} {input} > {output}"


# 6. run canu to assembly
rule canu_assembly:
    input:
        "iter_{n}.step_5.readHitSeq/union.fq"
    output:
        "iter_{n}.step_6.Canu/mito.contigs.fasta"
    params:
        canu_dir="iter_{n}.step_6.Canu"
    shell:
        "canu -p mito -d {params.canu_dir}  genomeSize={MITO_SIZE} -pacbio {input}"
