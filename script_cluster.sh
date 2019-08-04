# 设置 cluster 上可以使用的节点名称
CAPABLE_NODES="c1n1 c1n2 c1n3 c1n4 c1n5 c1n6 c2n2 c2n3 c2n4 c2n5 c2n6 rack1"

# 创建需要的文件夹
mkdir probes
mkdir results

# 这个终端函数用来在 cluster 上执行 alignment
# 因为之后的 MH catching 需要用到 clipped 序列信息 所以这里开启了 bwa mem 的 -Y 选项（全部使用 soft clipping）
alignForMTDCatching() {
    # 输入参数
    # 1 流程名称 用于给 cluster 作业以及输出文件起名
    # 2 reference fasta 文件
    # 3 fastq read 1 文件
    # 4 fastq read 2 文件
    bsub -J "$1-align" -m "$CAPABLE_NODES" -n 10 -R "span[ptile=10]" -q normal \
    "   bwa index $2
        bwa mem -t 8 -Y $2 $3 $4 | samtools view -@ 2 -F 2048 -Sb | samtools sort -@ 2 -m 10G -o $1.sort.bam
    "
}

# 这个函数将会执行四个步骤
# 1 查找序列中全部的 MH pairs 并生成 MTD 的特征序列
# 2 在 MH pairs 的关键位置统计 read depth
# 3 在 alignment 中查找具有 MTD 特征序列的 reads
# 4 使用 read depth 来标准化 reads count
countForSignatures() {
    # 输入参数
    # 1 流程名称 用于给 cluster 作业以及输出文件起名
    # 2 reference fasta 文件
    # 3 alignment map BAM 文件
    bsub -J "$1-depth" -m "$CAPABLE_NODES" -n 4 -R "span[ptile=4]" -q normal \
    "   find_mh probes/$1.mh < $2
        generate_signatures.awk $2 \`ls probes/$1.mh* | sort -t. -k1n,1\` > probes/signatures.$1.out
        samtools index -@ 4 $3
        sigtobed.awk probes/signatures.$1.out > probes/$1.pos.bed
        samtools depth -b probes/$1.pos.bed -d 0 $3 > mapped/$1.cov &
        samtools view $3 | catch_signatures.awk probes/signatures.$1.out > results/$1.sign.count.tsv &
        wait
        normalize_count.awk mapped/$1.cov results/$1.sign.count.tsv > results/$1.sign.norm.txt
    "
}


alignForMTDCatching lib1 ref/lib1.fa data/Library_1_clean_R1.fastq.gz data/Library_1_clean_R2.fastq.gz mapped/lib1.sort.bam
countForSignatures lib1 ref/lib1.fa mapped/lib1.sort.bam
