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
params.pdb = 'input.pdb'

process preparePDB {
    input:
    path input	

    output:
    path("${params.init}.gro")

    shell:
    '''
    gmx pdb2gmx -f ${input} -o ${params.init}.gro -water spce
    '''
}

process createNewBox {
    input:
    path input

    output:
    path("${params.init}_processed.gro")

    shell:
    '''
    gmx editconf -f ${params.init}.gro -o ${params.init}.gro -c -d 1.0 -bt cubic
    '''
}

process solvate {
    input:
    path input
    output:
    path("${params.init}_solvated.gro")
    path("topol.top")


    shell:
    '''
    gmx genbox -cp ${params.init}.gro -cs spc216.gro -o ${params.init}.gro -p topol.top
    '''
}

process ionize {
    input:
    file solvatedFile from solvate.solvatedFile
    output:
    file ionizedFile into '1D5R_solv_ions.gro'

    shell:
    '''
    gmx grompp -f ions.mdp -c ${solvatedFile} -p topol.top -o ions.tpr
    gmx genion -s ions.tpr -o ${ionizedFile} -p topol.top -pname NA -nname CL -nn 12
    '''
}


process minimization {

    input:
        path(topol_top_file)  // no `'` quotes here, just a variable name - the actual filename is handled in the workflow by the channels
        path(index_ndx_file)
    output:
        path("${params.mini_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.mini_prefix}.mdp -o ${params.mini_prefix}.tpr -c ${params.init}.gro -r ${params.init}.gro -p ${topol_top_file} -n ${index_ndx_file} -maxwarn -1
    gmx_d mdrun -v -deffnm ${params.mini_prefix}
    """
}

process equilibration {
    input:
        path("${params.mini_prefix}.gro")
        path(topol_top_file)
        path(index_ndx_file)
    output:
        path("${params.equi_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.equi_prefix}.mdp -o ${params.equi_prefix}.tpr -c ${params.mini_prefix}.gro -r ${params.init}.gro -p topol.top -n index.ndx
    gmx mdrun -v -deffnm ${params.equi_prefix}
    """
}

process production {
    input:
        path("${params.equi_prefix}.gro")
        path(topol.top)
        path(index.ndx)
    output:
        path("${params.prod_prefix}.tpr")
    script:
    """
    gmx grompp -f ${params.prod_prefix}.mdp -o ${params.prod_prefix}.tpr -c ${params.equi_prefix}.gro -p topol.top -n index.ndx
    gmx mdrun -v -deffnm ${params.prod_prefix}
    """
}

process figures {
    input:
        path("${params.mini_prefix}.edr")
        path("${params.equi_prefix}.edr")
        path("${params.prod_prefix}.edr")
    output:
        path("energy.xvg"),
        path("temp.xvg"),
        path("pressure.xvg")
    script:
    """
    echo 12 | gmx energy -f mini.edr -o energy.xvg
    echo 16 | gmx energy -f nvt.edr -o temp.xvg
    echo 17 | gmx energy -f npt.edr -o pressure.xvg
    """
}

workflow simulation {
    prepare_PDB_ch = Channel.fromPath(params.pdb)
    preparePDB(prepare_PDB_ch)
    
    
    

    ch1 = Channel.fromPath('topol.top')
    ch2 = Channel.fromPath('index.ndx')
    minimization(ch1, ch2)


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
