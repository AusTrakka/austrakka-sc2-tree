rule alignment_to_vcf:
    input:
        alignment=rules.nextclade.output.alignment,
    output:
        masked_vcf=temp('{group}.masked.vcf')
    params:
        mask = RESOURCES / "problematic_sites_sarsCov2.vcf",
        reference_sequence = RESOURCES / "MN908947.3.fna",
        reference_id = "MN908947.3",
    conda:
        ENVS / "usher.yaml"
    shell:
        """
        if ! grep -q "{params.reference_id}" "{input.alignment}"; then
            echo "Reference sequence {params.reference_id} not found in alignment {input.alignment}"
            echo "Adding reference sequence to alignment"
            cat {params.reference_sequence} >> {input.alignment}
        fi
        faToVcf -maskSites={params.mask} -ref={params.reference_id} {input.alignment} {output.masked_vcf}
        """

rule get_starting_tree:
    output:
        tree=temp('{group}.starting.nwk')
    params:
        starting=config["tree"].get("starting", False)
    conda:
        ENVS / "usher.yaml"
    shell:
        """
        if {params.starting}; then
            cp {params.starting} {output.tree}
        else 
            echo "()" > {output.tree}
        fi
        """

rule usher:
    input:
        vcf=rules.alignment_to_vcf.output.masked_vcf,
        starting_tree=rules.get_starting_tree.output.tree,
    output: 
        tree=temp('{group}.tree.pb'),
    params:
        reference = RESOURCES / "MN908947.3.fna"
    threads: 
        config["tree"].get("threads") if config["tree"].get("threads") else workflow.cores
    conda:
        ENVS / "usher.yaml"
    log:
       LOGS / "tree" / "usher.{group}.log"
    shell: 
        """
        usher \
          --threads {threads} \
          --collapse-tree \
          --tree {input.starting_tree} \
          --vcf {input.vcf} \
          -d {resources.tmpdir} \
          --save-mutation-annotated-tree {output.tree} 2>&1 | tee {log}
        """

rule matOptimize:
    input:
        tree=rules.usher.output.tree,
    output: 
        optimized_tree=temp('{group}.optimized.pb')
    threads: 
        config["tree"].get("threads") if config["tree"].get("threads") else workflow.cores
    conda:
        ENVS / "usher.yaml"
    log:
        LOGS / "tree" / "matOptimize.{group}.log"
    shell: 
        """
        matOptimize \
          --threads {threads} \
          --do-not-write-intermediate-files \
          -i {input.tree} \
          -o {output.optimized_tree} 2>&1 | tee {log}
        """

rule extract_tree:
    input:
        tree=rules.matOptimize.output.optimized_tree,
    output: 
        newick='{group}.nwk',
    threads: 
        config["tree"].get("threads") if config["tree"].get("threads") else workflow.cores
    conda:
        ENVS / "usher.yaml"
    log:
        LOGS / "tree" / "extract_tree.{group}.log"
    shell: 
        """
        matUtils extract -i {input.tree} -t {output.newick} 2>&1 | tee {log}
        """
