/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.seqware.workflows;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 *
 * @author pruzanov
 */
public class SampleFingerprintingWorkflow extends AbstractWorkflowDataModel {

    private String [] bam_files;
    private String [] genotypes;
    private String [] vcf_files;
    private String existingMatrix = "";
    private String studyName = "";
    private String watchersList = "";
    private String finalOutDir;
    private String dataDir;
    private String gatkVersion;
    private String tabixVersion;
    private String genomeFile;
    private String dbSNP_data;
    private String checkedSNPs;
    private String queue;
    
    //GATK parameters
    private String stand_call_conf = "50.0";
    private String stand_emit_conf = "10.0";
    private String dcov = "50";
    private int jChunkSize = 100; // Maximum allowed number of vcf files when jaccard_indexing step doesn't fork into multiple sub-jobs
    
    
    @Override
    public Map<String, SqwFile> setupFiles() {
     // Set up reference, bam and vcf files here     
     try {
       if (getProperty("bam_files") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "bam_files is not set, we need at least one bam file");
           return (null);
       } else {
           this.bam_files = getProperty("bam_files").split(",");
       }
       
       if (getProperty("genotypes") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "genotypes (vcf files) is not set, we need these to calculate similarity/generate genotypes");
           return (null);
       } else {
           this.genotypes = getProperty("genotypes").split(",");
       }
       
       if (getProperty("genome_file") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "genome_file (reference assembly fasta) is not set, we need it to generate a genotype");
           return (null);
       } else {
           this.genomeFile = getProperty("genome_file");
       }
       
       if (getProperty("dbsnp_file") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "dbsnp_file (reference dbsnp data) is not set, we need it to generate a genotype");
           return (null);
       } else {
           this.dbSNP_data = getProperty("dbsnp_file");
       }
       
       if (getProperty("checked_snps") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "checked_snps (checkpoints) is not set, we need it to generate a genotype");
           return (null);
       } else {
           this.checkedSNPs = getProperty("checked_snps");
       }
       
       if (getProperty("exisiting_matrix") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "exisiting_matrix not provided, will calclate jaccard indexes for all genotypes");
       } else {
           this.existingMatrix = getProperty("exisiting_matrix");
       }
       
       if (getProperty("watchers_list") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "watchers list not provided, will not send alerts");
       } else {
           this.watchersList = getProperty("watchers_list");
       }
       
       if (getProperty("stand_call_conf") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "stand_call_conf is not set, will use the default value of 50.0");
       } else {
           this.stand_call_conf = getProperty("stand_call_conf");
       }
       
       if (getProperty("stand_emit_conf") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "stand_emit_conf is not set, will use the default value of 10.0");
       } else {
           this.stand_emit_conf = getProperty("stand_emit_conf");
       }
       
       if (getProperty("dcov") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "dcov is not set, will use the default value of 50");
       } else {
           this.dcov = getProperty("dcov");
       } 
       
       if (getProperty("gatk_version") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "gatk_version is not set, we need it to call GATK correctly");
           return (null);
       } else {
           this.gatkVersion = getProperty("gatk_version");
       }
       
       if (getProperty("tabix_version") == null) {
           Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "tabix_version is not set, we need it to call tabix correctly");
           return (null);
       } else {
           this.tabixVersion = getProperty("tabix_version");
       }
           
       if (getProperty("queue") == null) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "Queue not set, most likely will run as default queue");
            this.queue = "";
        } else {
            this.queue = getProperty("queue");         
        }
       
        if (getProperty("study_name") == null) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "Study name isn't not set, will try to extract it from file names");
            this.studyName = "";
        } else {
            this.studyName = getProperty("study_name");         
        }
       
       // iterate over inputs
      for (int i = 0; i< this.bam_files.length; i++) {
        Log.stdout("CREATING FILE: bam_inputs_" + i);
        SqwFile file = this.createFile("bam_inputs_" + i);
        file.setSourcePath(this.bam_files[i]);
        file.setType("application/bam");
	file.setIsInput(true);
	Log.stdout("PROVISIONED path for bam_inputs_" + i + " is " + file.getProvisionedPath());
        
        //Using first file, try to guess the study name if it was not provided as an argument in .ini file
        if (i == 0 && this.studyName.isEmpty()) {
          if (bam_files[i].matches("SWID_\\d+_\\D+_\\d+")) {
             String tempName = bam_files[i].substring(bam_files[i].lastIndexOf("SWID_"));
             this.studyName = tempName.substring(0,tempName.indexOf("_") - 1);
           } else {
             this.studyName = "INDEF";
           }
        }
        
        SqwFile file_vcf = this.createFile("vcf_inputs_" + i);
        file_vcf.setType("text/vcf-4");
        String basename = this.bam_files[i].substring(this.bam_files[i].lastIndexOf("/")+1);
        String snpName  = basename.substring(0,basename.lastIndexOf(".bam")) + ".snps.raw.vcf";
        
        this.vcf_files = new String[this.bam_files.length];
        //If we don't have a vcf file we need to set it as an output
	if (null == this.genotypes[i] || this.genotypes[i].isEmpty() || this.genotypes[i].equals("NA")) {
          file_vcf.setSourcePath(this.dataDir + snpName);
          file_vcf.setIsOutput(true);
          file_vcf.setOutputPath(finalOutDir + snpName);
          file_vcf.setForceCopy(true);
          this.genotypes[i] = this.dataDir + snpName;
          // We will need to have a STUDYNAME_jaccard.matrix.csv output file setup here 
          // since we know by now that we will be generating new vcf files
          if (null==this.getFiles().get("jaccard_matrix")) {
           SqwFile file_matrix = this.createFile("jaccard_matrix");
           file_matrix.setType("text/plain");
           file_matrix.setSourcePath(this.dataDir + this.studyName + "_jaccard.matrix.csv");
           file_matrix.setIsOutput(true);
           file_matrix.setOutputPath(finalOutDir + this.studyName + "_jaccard.matrix.csv");
           file_matrix.setForceCopy(true);
          }
        } else {
          file_vcf.setSourcePath(this.genotypes[i]);
          file_vcf.setIsInput(true);
        }
        this.vcf_files[i] = file_vcf.getProvisionedPath();
      }
      
      // provision reference files
      SqwFile dbsnpFile = this.createFile("dbsnp_reference");
      dbsnpFile.setSourcePath(this.dbSNP_data);
      dbsnpFile.setType("application/vcf-4-gzip");
      dbsnpFile.setIsInput(true);
      
      //Set up the output (we have to set up only those which don't exist yet)
      SqwFile report_file = this.createFile("report_zip");
      String report_name = "sample_fingerprint";    
      report_file.setSourcePath(this.dataDir + report_name + ".report.zip");
      report_file.setType("application/zip-report-bundle");
      report_file.setIsOutput(true);
      report_file.setOutputPath(this.finalOutDir + report_name + ".report.zip");
      report_file.setForceCopy(true);
              
     } catch (Exception e) {
       Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, null, e);     
     }
     return this.getFiles();
    }
    
    @Override
    public void setupDirectory() {
        try {
          String outdir = null;
         if (getProperty("output_dir") != null && getProperty("output_prefix") != null) {
          outdir = getProperty("output_prefix") + getProperty("output_dir")
                             + "/seqware-" + this.getSeqware_version() + "_" + this.getName() + "_" + this.getVersion() + "/" + this.getRandom() + "/";
         }

         this.finalOutDir = outdir;
         this.addDirectory(getProperty("data_dir"));
         this.dataDir = getProperty("data_dir").endsWith("/") ? getProperty("data_dir") : getProperty("data_dir") + "/";
                 
       } catch (Exception e) {
         Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, null, e);
       }
        
    }
    
    
    @Override
    public void buildWorkflow() {
        try{
          //Need the decider to check the headers of the .bam files for the presence of RG and abscence of empty Fields (Field: )
          Workflow workflow = this.getWorkflow();
          int newVcfs = 0;
          List<Job> gatk_jobs = new ArrayList<Job>();
          
          for (int i = 0; i < this.bam_files.length; i++) {
           SqwFile vcf = this.getFiles().get("vcf_inputs_" + i);
           if (vcf.isInput()) {continue;} // Skip the files which were previously generated
           
           // First Step: Find which vcf files we need to generate
           if (null == this.genotypes[i] || this.genotypes[i].isEmpty() || this.genotypes[i].equals("NA")) {
             Job job_gatk = workflow.createBashJob("call_snps");
             job_gatk.setCommand("java -Xmx2g -jar " + getWorkflowBaseDir() + "/bin/GenomeAnalysisTK-" + this.gatkVersion + " GenomeAnalysisTK.jar "
                               + "-R " + this.genomeFile + " "
                               + "-T UnifiedGenotyper "
                               + "-I " + this.bam_files[i] + " "
                               + "--dbsnp " + this.getFiles().get("dbsnp_reference").getProvisionedPath() + " "
                               + "-o " + vcf.getSourcePath() + " "
                               + "-stand_call_conf " + this.stand_call_conf + " "
                               + "-stand_emit_conf " + this.stand_emit_conf + " "
                               + "-dcov " + this.dcov + " "
                               + "-L " + this.checkedSNPs);

            job_gatk.setMaxMemory(getProperty("gatk_memory"));
            if (!this.queue.isEmpty()) {
             job_gatk.setQueue(this.queue);
            }
            
            job_gatk.addFile(vcf);
            newVcfs++; // Used only to track the number of newly generated vcf files
            
            String basename = this.bam_files[i].substring(this.bam_files[i].lastIndexOf("/")+1);
            Job job_gatk2 = workflow.createBashJob("calculate_depth");
            job_gatk2.setCommand("java -Xmx3g -jar " + getWorkflowBaseDir() + "/bin/GenomeAnalysisTK-" + this.gatkVersion + " GenomeAnalysisTK.jar "
                            + "-R " + this.genomeFile + " "
                            + "-T DepthOfCoverage "
                            + "-I " + this.bam_files[i] + " "
                            + "-o " + this.dataDir + basename.substring(0,basename.lastIndexOf(".bam")) 
                            + "-L " + this.checkedSNPs);
            job_gatk2.setMaxMemory(getProperty("gatk_memory"));
            if (!this.queue.isEmpty()) {
              job_gatk2.setQueue(this.queue);
            }
            job_gatk2.addParent(job_gatk);
            
            //Additional step to create depth of coverage data - for individual fingerprint image generation
            Job job_fin = workflow.createBashJob("make_fingerprint_file");
            job_fin.setCommand(getWorkflowBaseDir() + "dependencies/create_fin.pl "
                          + "--refvcf="   + this.checkedSNPs + " "
                          + "--genotype=" + vcf.getSourcePath() + " " 
                          + "--coverage=" + this.dataDir + basename.substring(0,basename.lastIndexOf(".bam")) + ".sample_interval_summary "
                          + "--datadir="  + this.dataDir + " "
                          + "--basename=" + basename.substring(0,basename.lastIndexOf(".bam")));
            job_fin.setMaxMemory("3000");
            job_fin.addParent(job_gatk2);
            gatk_jobs.add(job_fin);         
            } 
           }
           
           // We don't need to continue if there are no new vcf files to generate
           if (newVcfs == 0) {
               Logger.getLogger(getClass().getName()).log(Level.INFO,"There are no new genotypes to generate, to avoid duplicate calculation the workflow will terminate");
               System.exit(0);
           }
           
           // At his point we need to stage, bgzip and tabix all vcf files we would need to create a jaccard matrix
           StringBuilder vcf_list    = new StringBuilder();
           vcf_list.append(vcf_files[0]);
           for (int i = 1;i < this.vcf_files.length; i++) {
               vcf_list.append(",").append(this.vcf_files[i]);
           }
           
           Job job_vcfprep = workflow.createBashJob("prepare_vcfs");
           job_vcfprep.setCommand("echo \"" + vcf_list.toString() + "\" "
                                + "perl -ne " + "\'{chomp;@vcfs=split(\",\");map{`"
                                + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip -c $_ > $_.gz && "
                                + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix -p vcf $_.gz`;} @vcfs}");
           job_vcfprep.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
             job_vcfprep.setQueue(this.queue);
            }
           
           //SET PARENTS FOR VCFPREP JOB
            for(Job parent: gatk_jobs) {
	       Log.stdout("Adding parent for tabix job");
	       job_vcfprep.addParent(parent);
	     }
           
           List<Job> upstream_jobs = new ArrayList<Job>();
           upstream_jobs.add(job_vcfprep);
           
           // We use simple reasoning here - if existingMatrix is unavailable, we will need to run a loooong job
           // of calculting jaccard indexes if the list of vcf files is long. So, we will break the list into smaller chunks
           // And use the results as exisitngMatrix in the next step Likely to run only once
           // (when running this workflow on a study for the very first time)
           if (this.existingMatrix.isEmpty() && this.vcf_files.length > this.jChunkSize) {
               // TODO:
               int currentChunkSize = this.jChunkSize;
               // Run x jobs on chunks of currentChunkSize files, recursively link to the next index calculator 
               int resultID = 1;
               while (this.vcf_files.length > currentChunkSize) {
               List<Job> jaccard_chunks = new ArrayList<Job>();
               StringBuilder vcf_chunk  = new StringBuilder();   // For list of vcf files in a chunk
               StringBuilder existingChunks = new StringBuilder(); // For the results from intermediate jchunk jobs
               
               
               int chunkMultiplier = 1;
               for (int vc = 0; vc < this.vcf_files.length; vc++) {
                 if (vc <= currentChunkSize * chunkMultiplier) {
                   // Accumulate jChunkSize files, set up job using previous results, if available
                   vcf_chunk.append(",").append(this.vcf_files[vc]); // jaccard_index script will ignore empty entries, comma at the start is Ok
                   continue;
                 }
                
                Job job_jchunk = workflow.createBashJob("make_matrix");
                String chunkName = "jaccard.chunk." + resultID + ".csv";
                resultID+=1;
                
                existingChunks.append(",").append(chunkName);
                job_jchunk.setCommand(getWorkflowBaseDir() + " jaccard_coeff.matrix.pl "
                            + "--list=" + vcf_chunk.toString() + " "
                            + "--studyname=" + this.studyName + " "
                            + "> " + this.dataDir + chunkName);
               if (!this.existingMatrix.isEmpty()) {
                  job_jchunk.getCommand().addArgument("--existing_matrix=" + this.existingMatrix);
               }
           
               job_jchunk.setMaxMemory("2000");
               if (!this.queue.isEmpty()) {
                job_jchunk.setQueue(this.queue);
               }
               
               for(Job parent: upstream_jobs) {
                job_jchunk.addParent(parent);
               } 
               chunkMultiplier += 1; 
               jaccard_chunks.add(job_jchunk);  
               vcf_chunk  = new StringBuilder();
               }
               // The result of the last job will be existingMatrix
               this.existingMatrix = existingChunks.toString();
               upstream_jobs.clear();
               upstream_jobs = jaccard_chunks;
               jaccard_chunks.clear();
                
               currentChunkSize+=this.jChunkSize;               
               }             
           }
           
           // Next job is jaccard_coeff.matrix script, need to create|update giant matrix of jaccard indexes
           // using vcftools, colors should not be assigned at his step. Will be final jaccard index job          
           SqwFile matrix = this.getFiles().get("jaccard_matrix");
           Job job_jaccard = workflow.createBashJob("make_matrix");
           job_jaccard.setCommand(getWorkflowBaseDir() + " jaccard_coeff.matrix.pl "
                                + "--list=" + vcf_list.toString() + " "
                                + "--studyname=" + this.studyName + " "
                                + "> " + matrix.getSourcePath());
           if (!this.existingMatrix.isEmpty()) {
             job_jaccard.getCommand().addArgument("--existing_matrix=" + this.existingMatrix);
           }
           
           job_jaccard.setMaxMemory("2000");
           if (!this.queue.isEmpty()) {
            job_jaccard.setQueue(this.queue);
           }
           job_jaccard.addFile(matrix);
           for(Job parent: upstream_jobs) {
              job_jaccard.addParent(parent);
           }
           
           // Images generated here: matrix slicing/color assignment, flagging suspicious files, wrapper for R script
           Job make_pics = workflow.createBashJob("make_report");
           make_pics.setCommand("perl " + getWorkflowBaseDir() + "dependencies/make_report.pl "
                              + "--matrix=" + matrix.getSourcePath() + " "
                              + "--datadir=" + this.dataDir + " "
                              + "--study=" + this.studyName);
           make_pics.setMaxMemory("2000");
           if (!this.queue.isEmpty()) {
            make_pics.setQueue(this.queue);
           }
           //TODO : make_report should write a message file if any heatmaps get flagged
           make_pics.addFile(matrix);
           make_pics.addParent(job_vcfprep);
        
           // Zip everything into a report bundle for provisioning
           Job zip_report = workflow.createBashJob("zip_everything");
           zip_report.setCommand("zip " + getFiles().get("report_zip").getSourcePath() + " "
                               + this.dataDir + "*.png "
                               + this.dataDir + matrix.getSourcePath() + " "
                               + this.dataDir + "*.html");
           
           zip_report.setMaxMemory("2000");
           zip_report.addParent(make_pics);
           zip_report.addFile(getFiles().get("report_zip"));
           
           // Next job is report-emailing script (check html for flagged samples, send email(s) to interested people)
           // Attach zipped report? Maybe NO
           //TODO : plotReporter should detect a file with a report on flagged files and send an email if needed
           if (!this.watchersList.isEmpty()) {
            Job job_alert = workflow.createBashJob("send_alert");
            job_alert.setCommand("perl" + getWorkflowBaseDir() + "dependencies/plotReporter.pl "
                              + "--emails=" + this.watchersList + " "
                              + "--report=" + this.dataDir + "index.html"
                              + "--zip=" + getFiles().get("report_zip").getSourcePath()
                              + "--finaldir=" + getFiles().get("report_zip").getOutputPath());

            job_alert.setMaxMemory("2000");
            job_alert.addParent(zip_report);
            if (!this.queue.isEmpty()) {
             job_alert.setQueue(this.queue);
            }
           }
           
          } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
          }
        
    }
    
}
