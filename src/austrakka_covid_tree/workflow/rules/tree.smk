rule alignment_to_vcf:
    input:
        alignment=rules.nextclade.output.alignment,
    output:
        masked_vcf=temp('{group}.masked.aln')
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
        newick='{group}.nwk',
        tree=temp('{group}.tree.pb'),
        optimized_tree=temp('{group}.optimized.pb')
    params:
        reference = RESOURCES / "MN908947.3.fna"
    threads: 
        config["tree"].get("threads") if config["tree"].get("threads") else workflow.cores
    conda:
        ENVS / "usher.yaml"
    log:
        LOGS / "tree.{group}.log"
    shell: 
        """
        usher \
          --collapse-tree \
          --tree {input.starting_tree} \
          --vcf {input.vcf} \
          -d {resources.tmpdir} \
          --save-mutation-annotated-tree {output.tree} > {log}

        matOptimize \
          --threads {threads} \
          --do-not-write-intermediate-files \
          -i {output.tree} \
          -o {output.optimized_tree} >> {log}
          
        matUtils extract -i optimized-tree.pb -t {output.newick} >> {log}
        """

