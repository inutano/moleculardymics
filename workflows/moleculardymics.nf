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


process preparePDB {
    input:
    path pdb_file

    output:
    path("${params.init}.gro"), emit: gro_input
    path 'topol.top', emit: topology

    shell:
    """
    echo 6 | gmx pdb2gmx -f ${pdb_file} -o ${params.init}.gro -water spce
    """
}

process createNewBox {
    input:
    path gro_file

    output:
    path("${params.init}_processed.gro"), emit: gro_processed

    shell:
    """
    gmx editconf -f ${gro_file} -o ${params.init}_processed.gro -c -d 1.0 -bt cubic
    """
}


process solvate {
    input:
    path gro_file
    path topol_file

    output:
    path("${params.init}_solvated.gro"), emit: gro_solvated
    path 'topol.top', emit: topology

    shell:
    """
    gmx solvate -cp ${gro_file} -cs spc216.gro -o ${params.init}_solvated.gro -p ${topol_file}
    """
}

process ionize {
    input:
    path gro_file
    path topol_file

    output:
    path("${params.init}_ionized.gro"), emit: gro_ionised
    path 'ions.tpr', emit: ions_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ions.mdp -c ${params.init}_solvated.gro -p ${topol_file} -o ions.tpr
    echo 13 | gmx genion -s ions.tpr -o ${params.init}_ionized.gro -p ${topol_file} -pname NA -nname CL -neutral
    """
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
    echo 12 | gmx energy -f ${params.mini_prefix}.edr -o energy.xvg
    echo 16 | gmx energy -f ${params.equi_prefix} -o temp.xvg
    echo 17 | gmx energy -f ${params.prod_prefix} -o pressure.xvg
    """
}

workflow {

    prepare_PDB_ch = Channel.fromPath("$baseDir/input.pdb")
    preparePDB(prepare_PDB_ch).set {PDB}
    preparePDB.out.topology.view()
    createNewBox(preparePDB.out.gro_input) 
    solvate(createNewBox.out.gro_processed, preparePDB.out.topology)
    ionize(solvate.out.gro_solvate, solavte.out.topology)
     
    }


// workflow simulation {
//    prepare_PDB_ch = Channel.fromPath(params.pdb)
 //   preparePDB(prepare_PDB_ch) | createNewBox | solvate | ionize | minimization | equilibration | production
    
//    createNewBox_ch = Channel.fromPath('topol.top')
//    ch2 = Channel.fromPath('index.ndx')
//   minimization(ch1, ch2)


 //   input:
 //       file('input.gro')
 //   output:
 //       file("energy.xvg"),
 //      file("temp.xvg"),
 //      file("pressure.xvg")
 //  minimization, equilibration, production, figures
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
