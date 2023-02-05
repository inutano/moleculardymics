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
params.run_prefix = 'md'



process preparePDB {
    input:
    path pdb_file

    output:
    path("${params.init}.gro"), emit: gro_input
    path 'topol.top', emit: topology
    path 'posre.itp', emit: posre_itp

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

process prepare_ionize {
    input:
    path ions_mdp
    path gro_file
    path topol_file

    output:
    path 'ions.tpr', emit: ions_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ${ions_mdp} -c ${params.init}_solvated.gro -p ${topol_file} -o ions.tpr
    """
}


process ionize {
    input:
    path ions_tplr
    path topol_file

    output:
    path("${params.init}_ionized.gro"), emit: gro_ionised
    path 'topol.top', emit: topology

    shell:
    """
    echo 13 | gmx genion -s ${ions_tpr} -o ${params.init}_ionized.gro -p ${topol_file} -pname NA -nname CL -neutral
    """
}

process prepare_minimise {
    input:
    path mini_mdp
    path gro_file
    path topol_file

    output:
    path("${params.mini_prefix}.tpr"), emit: mini_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ${mini_mdp} -c ${gro_file} -o ${params.mini_prefix}.tpr -p ${topol_file}
    """
}


process minimise {
    input:
    path mini_tpr
    path topol_file

    output:
    path("${params.mini_prefix}.gro"), emit: gro_mini
    path("${params.mini_prefix}.gro"), emit: edr_mini
    path("${params.mini_prefix}.gro"), emit: log_mini
    path("${params.mini_prefix}.gro"), emit: trr_mini
    path 'topol.top', emit: topology

    shell:
    """
    gmx mdrun -v -deffnm ${params.mini_prefix}
    """
}




process prepare_nvt {
    input:
    path nvt_mdp
    path gro_file
    path topol_file
    path posre_itp

    output:
    path("${params.equi_prefix}.tpr"), emit: nvt_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ${nvt_mdp} -c ${gro_file} -r ${gro_file} -o ${params.equi_prefix}.tpr -p ${topol_file}
    """
}


process nvt {
    input:
    path nvt_tpr
    path topol_file
    path posre_itp

    output:
    path("${params.equi_prefix}.gro"), emit: gro_nvt
    path("${params.equi_prefix}.gro"), emit: edr_nvt
    path("${params.equi_prefix}.gro"), emit: log_nvt
    path("${params.equi_prefix}.gro"), emit: trr_nvt
    path 'topol.top', emit: topology

    shell:
    """
    gmx mdrun -v -deffnm ${params.equi_prefix}
    """
}


process prepare_npt {
    input:
    path npt_mdp
    path gro_file
    path cpt_file
    path topol_file
    path posre_itp

    output:
    path("${params.prod_prefix}.tpr"), emit: npt_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ${npt_mdp} -c ${gro_file} -r ${gro_file} -t ${params.equi_prefix}.cpt -o ${params.prod_prefix}.tpr -p ${topol_file}
    """
}


process npt {
    input:
    path npt_tpr
    path topol_file

    output:
    path("${params.prod_prefix}.gro"), emit: gro_npt
    path("${params.prod_prefix}.gro"), emit: edr_npt
    path("${params.prod_prefix}.gro"), emit: log_npt
    path("${params.prod_prefix}.gro"), emit: trr_npt
    path 'topol.top', emit: topology

    shell:
    """
    gmx mdrun -v -deffnm ${params.prod_prefix}
    """
}


process prepare_md {
    input:
    path md_mdp
    path gro_file
    path cpt_file
    path topol_file
    path posre_itp

    output:
    path("${params.run_prefix}.tpr"), emit: md_tpr
    path 'topol.top', emit: topology

    shell:
    """
    gmx grompp -f ${md_mdp} -c ${gro_file} -r ${gro_file} -t ${params.prod_prefix}.cpt -o ${params.prod_prefix}.tpr -p ${topol_file}
    """
}


process md {
    input:
    path npt_tpr
    path topol_file

    output:
    path("${params.run_prefix}.gro"), emit: gro_md
    path("${params.run_prefix}.gro"), emit: edr_md
    path("${params.run_prefix}.gro"), emit: log_md
    path("${params.run_prefix}.gro"), emit: trr_md
    path 'topol.top', emit: topology

    shell:
    """
    gmx mdrun -v -deffnm ${params.prod_prefix}
    """
}


process figures {
    input:
        path("${params.mini_prefix}.edr")
        path("${params.equi_prefix}.edr")
        path("${params.prod_prefix}.edr")
        
    output:
        path("energy.xvg")
        path("temp.xvg")
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
    solvate(createNewBox.out.gro_processed, preparePDB.out.topology).set {solvate_ch}
   
    prepare_ions_ch = Channel.fromPath("$baseDir/ions.mdp")
    prepare_ionize(prepare_ions_ch, solvate.out.gro_solvated, solvate.out.topology)
    ionize(prepare_ionize.out.ions_tpr, prepare_ionize.out.topology)

    prepare_minimise_ch = Channel.fromPath("$baseDir/minim.mdp")
    prepare_minimise(prepare_minimise_ch, ionize.out.gro_ionised, ionize.out.topology)
    minimise(prepare_minimise.out.mini_tpr, prepare_minimise.out.topology)

    prepare_nvt_ch = Channel.fromPath("$baseDir/nvt.mdp")
    prepare_nvt(prepare_nvt_ch, minimise.out.gro_mini, minimise.out.topology, preparePDB.out.posre_itp)
    nvt(prepare_nvt.out.nvt_tpr, prepare_nvt.out.topology)

    prepare_npt_ch = Channel.fromPath("$baseDir/npt.mdp")
    prepare_npt(prepare_npt_ch, nvt.out.gro_nvt, nvt.out.cpt_file, nvt.out.topology, preparePDB.out.posre_itp)
    npt(prepare_npt.out.npt_tpr, prepare_npt.out.topology)

    prepare_md_ch = Channel.fromPath("$baseDir/md.mdp")
    prepare_md(prepare_md_ch, npt.out.gro_npt, npt.out.cpt_file, npt.out.topology, preparePDB.out.posre_itp)
    md(prepare_md.out.md_tpr, prepare_md.out.topology)

    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
