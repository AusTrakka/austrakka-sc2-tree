THREADS = int(config["tree"].get("threads")) if config["tree"].get("threads") else workflow.cores

rule alignment_to_vcf:
    """
    Converts the given alignment to vcf format and applies a mask to problematic sites.

    :input alignment:          The filtered alignment file produced by the :smk:ref:`filter` rule.

    :output masked_vcf:        The resulting masked vcf file.

    :param mask:               Path to the problematic sites vcf file in the RESOURCES directory that should be masked.
    :param reference_sequence: Path to the reference genome sequence (MN908947.3.fna) in the RESOURCES directory.
    :param reference_id:       The identifier for the reference genome sequence (MN908947.3).

    :conda:                    Path to the Conda environment file (usher.yaml) in the ENVS directory.

    .. note::
        This rule first checks whether the reference sequence ID is found in the given alignment. If not, it appends
        the reference sequence to the alignment file. Then, it converts the alignment to vcf format using the `faToVcf` 
        command, applying a mask to problematic sites. The problematic sites are defined in the provided mask vcf file.
    """
    input:
        alignment=rules.filter_nextclade.output.alignment_filtered,
    output:
        masked_vcf=temp('{outdir}/{name}.masked.vcf')
    params:
        mask = RESOURCES / "problematic_sites_sarsCov2.vcf",
        reference_sequence = RESOURCES / "MN908947.3.fna",
        reference_id = "MN908947.3",
        include_reference = 1 if config["tree"].get("include_reference", True) else 0
    conda:
        ENVS / "usher.yaml"
    shell:
        """
        if ! grep -q "^>{params.reference_id}$" "{input.alignment}"; then
            echo "Reference sequence {params.reference_id} not found in alignment {input.alignment}"
            echo "Adding reference sequence to alignment"
            cat {params.reference_sequence} >> {input.alignment}
        fi
        if [ ! -f {params.include_reference} ]; then
            cat {params.reference_sequence} >> {input.alignment}
        fi
        faToVcf -maskSites={params.mask} -ref={params.reference_id} {input.alignment} {output.masked_vcf}
        """

rule get_starting_tree:
    """
    Retrieves or creates a starting tree for the phylogenetic tree construction process.

    :output tree:              The starting tree in Newick format.

    :param starting:           The optional path to the starting tree. This path is fetched from the "starting" field
                               in the "tree" section of the configuration. If not specified, an empty tree is used.

    :conda:                    Path to the Conda environment file (usher.yaml) in the ENVS directory.

    .. note::
        If a starting tree is specified in the configuration, this rule simply copies the specified file to the output 
        location. If no starting tree is specified, it creates an empty tree (represented as "()") and writes it to the 
        output location.
    """
    output:
        tree=temp('{outdir}/{name}.starting.nwk')
    params:
        starting=config.get('starting_tree') if config.get('starting_tree', False) else 0
    conda:
        ENVS / "usher.yaml"
    shell:
        """
        if [ {params.starting} = 0 ]; then
            echo "()" > {output.tree}
        else
            cp {params.starting} {output.tree}
        fi
        """

rule usher:
    """
    Constructs a phylogenetic tree using UShER based on given vcf and starting tree.
    
    .. note::
        The UShER tool is used to quickly place new sequences into an existing phylogenetic tree while maintaining 
        the overall tree topology. In this rule, we are running UShER with several options including multi-threading,
        tree collapsing, and mutation-annotated tree saving. The standard error and output are logged in the specified 
        log file.

    :input vcf:               The masked vcf file produced by the :smk:ref:`alignment_to_vcf` rule.
    :input starting_tree:     The starting tree file produced by the :smk:ref:`get_starting_tree` rule.

    :output tree:             The constructed phylogenetic tree in protobuf format.

    :param reference:         Path to the reference genome file (MN908947.3.fna) in the RESOURCES directory.

    :threads:                 Number of threads to use for this rule, which is fetched from the "tree" section of the 
                            configuration. If not specified, the total number of available cores in the workflow 
                            is used.

    :conda:                   Path to the Conda environment file (usher.yaml) in the ENVS directory.

    :log:                     Log file path in the "tree" subdirectory of the LOGS directory. The filename is 
                            structured as "usher.{name}.log".

    """
    input:
        vcf=rules.alignment_to_vcf.output.masked_vcf,
        starting_tree=rules.get_starting_tree.output.tree,
    output: 
        tree='{outdir}/{name}.pb' if config["mat"] else temp('{outdir}/{name}.pb')
    params:
        batch_size_per_process=config["tree"].get("batch_size_per_process", 10),
        optimization_radius=config["tree"].get("optimization_radius", 0),
    threads: 
        THREADS
    conda:
        ENVS / "usher.yaml"
    log:
       LOGS / "tree" / "usher.{name}.log"
    shell: 
        """
        usher-sampled \
          --threads {threads} \
          --optimization_radius {params.optimization_radius} \
          --batch_size_per_process {params.batch_size_per_process} \
          --sort-before-placement-3 \
          --tree {input.starting_tree} \
          --vcf {input.vcf} \
          -d {resources.tmpdir} \
          --save-mutation-annotated-tree {output.tree} 2>&1 | tee {log}
        """

rule matOptimize:
    """
    Optimizes a Mutation Annotated Tree (MAT) using the matOptimize command from UShER.

    :input tree:               The phylogenetic tree in protobuf format produced by the :smk:ref:`usher` rule.

    :output optimized_tree:    The optimized MAT in protobuf format.

    :threads:                  Number of threads to use for this rule, which is fetched from the "tree" section of the 
                               configuration. If not specified, the total number of available cores in the workflow 
                               is used.

    :conda:                    Path to the Conda environment file (usher.yaml) in the ENVS directory.

    :log:                      Log file path in the "tree" subdirectory of the LOGS directory. The filename is 
                               structured as "matOptimize.{name}.log".

    .. note::
        The matOptimize command is used to optimize an input MAT by improving the parsimony score. The tool does not 
        write intermediate files during optimization (--do-not-write-intermediate-files option). The standard error and 
        output from the command are logged in the specified log file.
    """
    input:
        tree=rules.usher.output.tree,
    output: 
        optimized_tree=temp('{outdir}/{name}.optimized.pb')
    threads: 
        THREADS
    conda:
        ENVS / "usher.yaml"
    log:
        LOGS / "tree" / "matOptimize.{name}.log"
    shell: 
        """
        matOptimize \
          --threads {threads} \
          --do-not-write-intermediate-files \
          -i {input.tree} \
          -o {output.optimized_tree} 2>&1 | tee {log}
        """

rule extract_tree:
    """
    Extracts a tree from the optimized Mutation Annotated Tree (MAT) in Newick format using the matUtils command from UShER.

    :input tree:               The optimized MAT in protobuf format produced by the :smk:ref:`matOptimize` rule.

    :output newick:            The extracted tree in Newick format.

    :threads:                  Number of threads to use for this rule, which is fetched from the "tree" section of the 
                               configuration. If not specified, the total number of available cores in the workflow 
                               is used.

    :conda:                    Path to the Conda environment file (usher.yaml) in the ENVS directory.

    :log:                      Log file path in the "tree" subdirectory of the LOGS directory. The filename is 
                               structured as "extract_tree.{name}.log".

    .. note::
        The matUtils extract command is used to extract a tree from the input MAT and write it in Newick format. The 
        standard error and output from the command are logged in the specified log file.
    """
    input:
        tree=rules.matOptimize.output.optimized_tree,
    output: 
        newick='{outdir}/{name}.nwk',
    threads: 
        THREADS
    conda:
        ENVS / "usher.yaml"
    log:
        LOGS / "tree" / "extract_tree.{name}.log"
    shell: 
        """
        matUtils extract -i {input.tree} -t {output.newick} 2>&1 | tee {log}
        """
