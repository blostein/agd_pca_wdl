version 1.0

#IMPORTS
## According to this: https://cromwell.readthedocs.io/en/stable/Imports/ we can import raw from github
## so we can make use of the already written WDLs provided by WARP/VUMC Biostatistics

import "https://raw.githubusercontent.com/shengqh/warp/develop/tasks/vumc_biostatistics/GcpUtils.wdl" as http_GcpUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/genotype/Utils.wdl" as http_GenotypeUtils
import "https://raw.githubusercontent.com/shengqh/warp/develop/pipelines/vumc_biostatistics/agd/AgdUtils.wdl" as http_AgdUtils


workflow VUMCGenotypePCA {
  input {
    Array[File] source_pgen_files
    Array[File] source_pvar_files
    Array[File] source_psam_files

    Array[String] chromosomes

    String target_prefix

    File variants_extract_file
    File id_map_file
    File pca_loadings_file
    File pca_af_file

    String? project_id
    String? target_gcp_folder
  }

  scatter (idx in range(length(chromosomes))) {
    String chromosome = chromosomes[idx]
    File pgen_file = source_pgen_files[idx]
    File pvar_file = source_pvar_files[idx]
    File psam_file = source_psam_files[idx]
    String replaced_sample_name = "~{chromosome}.psam"

    #I think I need this to get the IDs correctly as GRIDS

    call http_AgdUtils.ReplaceICAIdWithGrid as ReplaceICAIdWithGrid {
      input:
        input_psam = psam_file,
        id_map_file = id_map_file,
        output_psam = replaced_sample_name
    }

    call ExtractVariants as ExtractVariants{
      input:
        pgen_file = pgen_file,
        pvar_file = pvar_file,
        psam_file = ReplaceICAIdWithGrid.output_psam,
        chromosome = chromosome,
        variants_extract_file = variants_extract_file
    }
  }

  call http_GenotypeUtils.MergePgenFiles as MergePgenFiles{
    input:
      pgen_files = ExtractVariants.output_pgen_file,
      pvar_files = ExtractVariants.output_pvar_file,
      psam_files = ExtractVariants.output_psam_file,
      target_prefix = target_prefix
  }

  call ProjectPCA{
    input: 
      pgen_file = MergePgenFiles.output_pgen_file, 
      pvar_file = MergePgenFiles.output_pvar_file,
      psam_file = MergePgenFiles.output_psam_file,  
      PCA_loadings = pca_loadings_file,
      PCA_AF = pca_af_file,
      OUTNAME = target_prefix
  }

  if(defined(target_gcp_folder)){
    call http_GcpUtils.MoveOrCopyOneFile as CopyFile_one {
      input:
        source_file = ProjectPCA.output_pca_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pca_file = select_first([CopyFile_one.output_file, ProjectPCA.output_pca_file])
  }
}

### TASK DEFINITIONS 

task ExtractVariants{
  input {
    File pgen_file
    File pvar_file
    File psam_file 

    String chromosome

    File variants_extract_file

    Int memory_gb

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  Int disk_size = ceil(size([pgen_file, psam_file, pvar_file], "GB")  * 2) + 20

  String new_pgen = chromosome + ".pgen"
  String new_pvar = chromosome + ".pvar"
  String new_psam = chromosome + ".psam"
  String intermediate_pgen = chromosome + "_varids.pgen"
  String intermediate_pvar = chromosome + "_varids.pvar"
  String intermediate_psam = chromosome + "_varids.psam"

  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      --snps-only \
      --set-all-var-ids @:#:\$r:\$a \
      --new-id-max-allele-len 1000 \
      --make-pgen \
      --out ~{chromosome}_varids
    
    plink2 \
      --pgen ~{intermediate_pgen} \
      --pvar ~{intermediate_pvar} \
      --psam ~{intermediate_psam} \
      --extract ~{variants_extract_file} \
      --make-pgen \
      --out ~{chromosome}
  }

  runtime {
    docker: docker
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }

  output {
    File output_pgen_file = new_pgen
    File output_pvar_file = new_pvar
    File output_psam_file = new_psam
  }

}

task ProjectPCA{
  input{
    File pgen_file
    File pvar_file
    File psam_file
    File PCA_loadings
    File PCA_AF
    String OUTNAME

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

  }

  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB")  * 2) + 20

  String pca_file = OUTNAME + ".genotype.pca.sscore"

  command {
    plink2 --pgen ~{pgen_file} --pvar ~{pvar_file} --psam ~{psam_file} --score ~{PCA_loadings} \
    variance-standardize \
    cols=-scoreavgs,+scoresums \
    list-variants \
    header-read \
    --score-col-nums 3-12 \
    --read-freq ~{PCA_AF} \
    --out ~{OUTNAME}

    cp ~{OUTNAME}.sscore ~{pca_file}
    }

  runtime {
    docker: docker
    preemptible: 1
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    memory: memory_gb + " GiB"
  }
  
  output {
    File output_pca_file = "~{pca_file}"
  }
}


