#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.init = 'input'
params.mini_prefix = 'mini'
params.equi_prefix = 'nvt'
params.prod_prefix = 'npt'


process preparePDB {
    conda: 'env.yaml'
    
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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

    input:
    path ions_tpr
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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

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
    conda: 'env.yaml'

    input:
    path nvt_tpr
    path topol_file

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


    }
