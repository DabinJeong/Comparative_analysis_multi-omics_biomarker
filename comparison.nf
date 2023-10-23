#! /usr/bin/env nextflow
/*
 * Enable dsl2
 */
nextflow.enable.dsl = 2

params.cpus = "30"
params.mem = "10" //GB

process TrainTestSplit{
        publishDir "${params.publish_dir}_${params.iter}_${params.kFold}/${task.process.replaceAll(':', '_')}", mode: "copy"
        
        input:
               path metabolome
               path methylome
               path proteome
               path transcriptome
               path clinical 
        output:
               publishDir 
               tuple file("${params.kFold}_${params.label}_samples_train.txt"),file("${params.kFold}_${params.label}_samples_test.txt"),file("${params.kFold}_Merged_clinical_${params.label}_train.tsv") 
        script:
        """
        cp $baseDir/result_20220912/ourBiomarker_GCN_out_${params.iter}/test_${params.kFold}/TrainTestSplit/* ./
        """
}

process DIABLO{
        publishDir "${params.publish_dir}_${params.iter}_${params.kFold}/${task.process.replaceAll(':', '_')}", mode: "copy"
        
        cpus params.cpus
        memory "${params.mem} GB" 
        
        input:
                path methylome
                path proteome
                path transcriptome
                path clinical
                tuple file(train_samples),file(test_samples),file(clinical_train) 
        output:
                publishDir
                file("DIABLO_performance.txt")
        
        script:
        """
        Rscript $baseDir/modules_comparativeAnalysis/DIABLO/diablo.R ${params.label} ${methylome} ${proteome} ${transcriptome} ${clinical} ${train_samples} ${test_samples}
        """
}

process MOGONET{
        publishDir "${params.publish_dir}_${params.iter}_${params.kFold}/${task.process.replaceAll(':', '_')}", mode: "copy"

        cpus params.cpus
        memory "${params.mem} GB"

        input:
                path methylome
                path proteome
                path transcriptome
                path clinical
                tuple file(train_samples),file(test_samples),file(clinical_train)
        output:
                publishDir
                file("MOGONET_performance.txt")

        script:
        """
        python $baseDir/modules_comparativeAnalysis/MOGONET/preprocessing.py --label ${params.label} --m ${methylome} --p ${proteome} --t ${transcriptome} --clin ${clinical} --train_samples ${train_samples} --test_samples ${test_samples}  
        python $baseDir/modules_comparativeAnalysis/MOGONET/main_mogonet.py
        """
}


process MOFA{
        publishDir "${params.publish_dir}_${params.iter}_${params.kFold}/${task.process.replaceAll(':', '_')}", mode: "copy"

        cpus params.cpus
        memory "${params.mem} GB"

        input:
                path methylome
                path proteome
                path transcriptome
                path clinical
                tuple file(train_samples),file(test_samples),file(clinical_train)
        output:
                publishDir
                file("MOFA_performance.txt")

        script:
        """
        Rscript $baseDir/modules_comparativeAnalysis/MOFA/MOFA.R ${params.label} ${methylome} ${proteome} ${transcriptome} ${clinical} ${train_samples} ${test_samples}
        """
}


process iDRW{
        publishDir "${params.publish_dir}_${params.iter}_${params.kFold}/${task.process.replaceAll(':', '_')}", mode: "copy"
        
        cpus params.cpus
        memory "${params.mem} GB" 
        
        input:
                path methylome
                path proteome
                path transcriptome
                path clinical
                tuple file(train_samples),file(test_samples),file(clinical_train) 
        output:
                publishDir
                file("iDRW_performance.txt")
        
        script:
        """
        Rscript $baseDir/modules_comparativeAnalysis/iDRW/iDRW.R ${params.label} ${methylome} ${proteome} ${transcriptome} ${clinical} ${train_samples} ${test_samples}
        """
}

/*
 * run pipeline
 */

workflow{
        TrainTestSplit(file(params.metabolome),
                       file(params.methylome),
                       file(params.proteome),
                       file(params.transcriptome),
                       file(params.clinical))    
        DIABLO(file(params.methylome),
               file(params.proteome),
               file(params.transcriptome),
               file(params.clinical), 
               TrainTestSplit.out)
        MOGONET(file(params.methylome),
                file(params.proteome),
                file(params.transcriptome),
                file(params.clinical),
                TrainTestSplit.out)
        iDRW(file(params.methylome),
               file(params.proteome),
               file(params.transcriptome),
               file(params.clinical), 
               TrainTestSplit.out)
        MOFA(file(params.methylome),
               file(params.proteome),
               file(params.transcriptome),
               file(params.clinical), 
               TrainTestSplit.out)
}
