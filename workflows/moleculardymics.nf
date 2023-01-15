/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMoleculardymics.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.init = 'input'
params.mini_prefix = 'mini'
params.equi_prefix = 'nvt'
params.prod_prefix = 'npt'

process minimization {
    input:
        file('topol.top'),
        file('index.ndx')
    output:
        file("${params.mini_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.mini_prefix}.mdp -o ${params.mini_prefix}.tpr -c ${params.init}.gro -r ${params.init}.gro -p topol.top -n index.ndx -maxwarn -1
    gmx_d mdrun -v -deffnm ${params.mini_prefix}
    """
}

process equilibration {
    input:
        file("${params.mini_prefix}.gro"),
        file('topol.top'),
        file('index.ndx')
    output:
        file("${params.equi_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.equi_prefix}.mdp -o ${params.equi_prefix}.tpr -c ${params.mini_prefix}.gro -r ${params.init}.gro -p topol.top -n index.ndx
    gmx mdrun -v -deffnm ${params.equi_prefix}
    """
}

process production {
    input:
        file("${params.equi_prefix}.gro"),
        file('topol.top'),
        file('index.ndx')
    output:
        file("${params.prod_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.prod_prefix}.mdp -o ${params.prod_prefix}.tpr -c ${params.equi_prefix}.gro -p topol.top -n index.ndx
    gmx mdrun -v -deffnm ${params.prod_prefix}
    """
}

process figures {
    input:
        file("${params.mini_prefix}.edr"),
        file("${params.equi_prefix}.edr"),
        file("${params.prod_prefix}.edr")
    output:
        file("energy.xvg"),
        file("temp.xvg"),
        file("pressure.xvg")
    script:
    """
    echo 12 | gmx energy -f mini.edr -o energy.xvg
    echo 16 | gmx energy -f nvt.edr -o temp.xvg
    echo 17 | gmx energy -f npt.edr -o pressure.xvg
    """
}

workflow simulation {
    input:
        file('input.gro')
    output:
        file("energy.xvg"),
        file("temp.xvg"),
        file("pressure.xvg")
    minimization, equilibration, production, figures
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
