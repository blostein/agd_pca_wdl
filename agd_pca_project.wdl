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
        outputName = chromosome,
        extractFile = variants_extract_file
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
      pfile_prefix = MergePgenFiles.output_pfile_prefix,
      PCA_loadings = pca_loadings_file,
      PCA_AF = pca_af_file,
      OUTNAME = target_prefix
  }

  if(defined(target_gcp_folder)){
    call http_GcpUtils.MoveOrCopyThreeFiles as CopyFile {
      input:
        source_file1 = MergePgenFiles.output_pgen_file,
        source_file2 = MergePgenFiles.output_pvar_file,
        source_file3 = MergePgenFiles.output_psam_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  if(defined(target_gcp_folder)){
    call GcpUtils.MoveOrCopyOneFile as CopyFile {
      input:
        source_file = ProjectPCA.output_pca_file,
        is_move_file = false,
        project_id = project_id,
        target_gcp_folder = select_first([target_gcp_folder])
    }
  }

  output {
    File output_pca_file = select_first([CopyFile.output_file, PlinkPCA.output_pca_file])
  }
}

### TASK DEFINITIONS 

task ExtractVariants{
  input {
    File pgen_file
    File pvar_file
    File psam_file 

    String outputName
    File variants_extract_file

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"
  }

  command {
    plink2 \
      --pgen ~{pgen_file} \
      --pvar ~{pvar_file} \
      --psam ~{psam_file} \
      --extract ~{variants_extract_file} \
      --make-pfile \
      --out ~{outputName}
  }

  runtime {
    docker: docker
  }

  output {
    File output_pgen_file = "~{outputName}.pgen"
    File output_pvar_file = "~{outputName}.pvar"
    File output_psam_file = "~{outputName}.psam"
  }

}

task ProjectPCA{
  input{
    String pfile_prefix
    File PCA_loadings
    File PCA_AF
    String OUTNAME

    Int memory_gb = 20
    Int cpu = 8

    String docker = "hkim298/plink_1.9_2.0:20230116_20230707"

  }

  String pgen_file = pfile_prefix + '.pgen'
  String pvar_file = pfile_prefix + '.pvar'
  String psam_file = pfile_prefix + '.psam'

  Int disk_size = ceil(size([pgen_file, pvar_file, psam_file], "GB")  * 2) + 20

  String pca_file = OUTNAME + ".genotype.pca.sscore"

  command {
    plink2 --pfile ~{pfile_prefix} --score ~{PCA_loadings} \
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


