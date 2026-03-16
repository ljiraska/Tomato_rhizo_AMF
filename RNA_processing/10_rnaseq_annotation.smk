### SORTME RNA
rule sortmerna:
    input:
        r1="/home/amhorst/2021_tomato_rhizo/RNA/reads/hostrm/{sample}_R1_rmhost.fq.gz",
        r2="/home/amhorst/2021_tomato_rhizo/RNA/reads/hostrm/{sample}_R2_rmhost.fq.gz",
        #cleanup="logs/cleanup/{sample}.done"
    output:
        #check_sortmerna="Checks/SortmeRNA/{sample}_sortmeRNA.done",
        check_repair = "Checks/SortmeRNA/{sample}_repair.done",
        r1_cleaned="outputs/sortmeRNA/{sample}/{sample}_unaligned_R1.fq.gz",
        r2_cleaned="outputs/sortmeRNA/{sample}/{sample}_unaligned_R2.fq.gz"
    resources:
        mem=24000,
        time="48:00:00",
        partition="med2"
    threads: 36
    params:
        ref="/group/jbemersogrp/users/lucie/DATABASES/sortmerna/rRNA_databases_v4",
        output_dir="/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/sortmeRNA",
        tag="{sample}"
    shell:
        """
        echo 'Running sortmeRNA'

        mkdir -p {params.output_dir}/{params.tag}

        pixi run sortmerna \
            --ref {params.ref}/smr_v4.3_default_db.fasta \
            --reads {input.r1} --reads {input.r2} \
            --workdir {params.output_dir}/{params.tag} \
            --aligned {params.output_dir}/{params.tag}/aligned \
            --other {params.output_dir}/{params.tag}/unaligned --fastx && \
            touch {output.check_sortmerna}

        echo 'Splitting interleaved unaligned reads with repair.sh'

        pixi run repair.sh \
            in={params.output_dir}/{params.tag}/unaligned.fq.gz \
            out1={output.r1_cleaned} \
            out2={output.r2_cleaned} \
            overwrite=true && \
            touch {output.check_repair}
        """

### Salmon indexing and quanitfication

rule salmon_index:
    input:
        assembly = "/home/amhorst/2021_tomato_rhizo/RNA/co_assemblies/first_co_ass/final_contigs/230125_all_coas_nucl.fa"
    output:
        check_salmon_index = "Checks/salmon_index/salmon_index.done"
    resources:
        mem = 24000,
        time = "24:00:00",
        partition = "med2"
    threads: 32
    params:
        output_dir = "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/salmon_index",
        tag = "salmon_index"
    shell:
        """
        echo 'Running Salmon index'

        mkdir -p {params.output_dir}

        pixi run salmon index -t {input.assembly} -i {params.output_dir} --threads {threads} && \
            touch {output.check_salmon_index}
        """

rule salmon_quantification:
    input:
        salmon_index_done = "Checks/salmon_index/salmon_index.done",
        check_repair = "Checks/SortmeRNA/{sample}_repair.done",
    output:
        check_salmon_quant = "Checks/salmon_index/salmon_quant_{sample}.done",
    resources:
        mem = 24000,
        time = "24:00:00",
        partition = "med2"
    threads: 32
    params:
        output_dir = "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/salmon_quant",
        tag="{sample}",
        salmon_index = "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/salmon_index",
        #salmon_index=lambda wildcards: "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/salmon_index",
        r1_clean="outputs/sortmeRNA/{sample}/{sample}_unaligned_R1.fq.gz",
        r2_clean="outputs/sortmeRNA/{sample}/{sample}_unaligned_R2.fq.gz"
    shell:
        """
        mkdir -p {params.output_dir}
        pixi run salmon quant -i {params.salmon_index} -l A \
                    -1 {params.r1_clean} -2 {params.r2_clean} \
                    -p {threads} \
                    -o {params.output_dir}/{params.tag} && \
            touch {output.check_salmon_quant}
        """

rule salmon_quantmerge:
    input:
        expand("outputs/salmon_quant/{sample}/quant.sf", sample=samples)
    output:
        tpm = "outputs/salmon_quant/merged_quant_TPM.tsv",
        numreads = "outputs/salmon_quant/merged_quant_numreads.tsv",
        check_salmon_merge_tpm = "Checks/salmon_merge/salmon_merge_tpm.done",
        check_salmon_merge_counts = "Checks/salmon_merge/salmon_merge_counts.done",
    resources:
        mem = 8000,
        time = "2:00:00",
        partition = "med2",
    params:
        tag="salmon_quantmerge",
    shell:
        """
        quant_dirs=$(ls -d outputs/salmon_quant/* | tr '\n' ',' | sed 's/,$//')
        echo "$quant_dirs"

        sample_names=$(ls -d outputs/salmon_quant/* | xargs -n1 basename | tr '\n' ',' | sed 's/,$//')
        echo "$sample_names"

        pixi run salmon quantmerge \
            --quants $quant_dirs \
            --names $sample_names \
            -c TPM \
            -o {output.tpm} && \
        touch {output.check_salmon_merge_tpm}

        pixi run salmon quantmerge \
            --quants $quant_dirs \
            --names $sample_names \
            -c numreads \
            -o {output.numreads} && \
        touch {output.check_salmon_merge_counts}
        """

### Transdecoder

rule transdecoder_predict:
    input:
        assembly = "/home/amhorst/2021_tomato_rhizo/RNA/co_assemblies/first_co_ass/final_contigs/230125_all_coas_nucl.fa"
    output:
        check_LongOrfs = "Checks/transdecoder/check_LongOrfs.done",
        check_Predict = "Checks/transdecoder/check_Predict.done"
    resources:
        mem = 24000,
        time = "36:00:00",
        partition = "med2"
    threads: 32
    params:
        output_dir = "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/transdecoder",
        tag = "transdecoder"
    shell:
        """
        mkdir -p {params.output_dir}

        #cp {input} transdecoder/assembly.fa

        pixi run TransDecoder.LongOrfs -t {input.assembly} --output_dir {params.output_dir} && \
            touch {output.check_LongOrfs}

        pixi run TransDecoder.Predict -t {input.assembly} --output_dir {params.output_dir} && \
            touch {output.check_Predict}
        """

### EggNOG mapper

rule eggnog_mapper:
    input:
        check_LongOrfs = "Checks/transdecoder/check_LongOrfs.done",
        check_Predict = "Checks/transdecoder/check_Predict.done"
    output:
        check_eggnog_mapper = "Checks/eggnog_mapper/check_eggnog_mapper.done"
    resources:
        mem = 24000,
        time = "36:00:00",
        partition = "med2"
    threads: 32
    params:
        output_dir = "/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/eggnog_mapper",
        longest_orfs ="/home/ljiraska/emersongrp_lucie/Projects/AMF/outputs/transdecoder/230125_all_coas_nucl.fa.transdecoder_dir/longest_orfs.pep",
        tag = "eggnog_mapper"
    shell:
        """
        mkdir -p {params.output_dir}
        pixi run emapper.py \
            --data_dir ../../DATABASES/eggnog \
            -i {params.longest_orfs} \
            -o {params.tag} \
            --output_dir {params.output_dir} \
            --itype proteins \
            --excel \
            --cpu {threads} && \
        touch {output.check_eggnog_mapper}
        """

