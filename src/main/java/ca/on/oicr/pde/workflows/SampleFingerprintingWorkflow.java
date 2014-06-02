/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 *
 * @author pruzanov
 */
public class SampleFingerprintingWorkflow extends OicrWorkflow {

    private String[] bam_files;
    private String[] genotypes;
    private String[] vcf_files;
    private String[] GATK_dirs;
    private String existingMatrix = "";
    private String studyName = "";
    private String watchersList = "";
    private String dataDir;
    private String tempDir = "tempfiles/";
    //Additional one for GATK:
    private String tmpDir = "temp";
    private String finDir = "finfiles/";
    private String gatkVersion;
    private String tabixVersion;
    private String vcftoolsVersion;
    private String samtoolsVersion;
    private String genomeFile;
    private String checkedSNPs;
    private String check_points;
    private String queue;
    private String reportName = "sample_fingerprint";
    private String gatk_java;
    //GATK parameters
    private String stand_call_conf = "50.0";
    private String stand_emit_conf = "10.0";
    private String dcov = "200";
    private final int jChunkSize = 50; // Maximum allowed number of vcf files when jaccard_indexing step doesn't fork into multiple sub-jobs
    private final int batchCount = 100; // Use for job batching, this many jobs 
    private boolean manualOutput;
    private boolean provisionVcfs;

    @Override
    public Map<String, SqwFile> setupFiles() {
        // Set up reference, bam and vcf files here     
        try {
            if (getProperty("input_files") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "input_files is not set, we need at least one bam file");
                return (null);
            } else {
                this.bam_files = getProperty("input_files").split(",");
            }

            if (getProperty("genotypes") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "genotypes (vcf files) is not set, we need these to calculate similarity/generate genotypes");
                return (null);
            } else {
                this.genotypes = getProperty("genotypes").split(",");
                if (this.genotypes.length <= 1) {
                    this.genotypes = null;
                }
            }

            if (getProperty("gatk_java") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "gatk_java is not set, we need it to run GATK since it requires 1.7 java and we cannot rely on the defaul");
                return (null);
            } else {
                this.gatk_java = getProperty("gatk_java");
            }

            if (getProperty("genome_file") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "genome_file (reference assembly fasta) is not set, we need it to generate a genotype");
                return (null);
            } else {
                this.genomeFile = getProperty("genome_file");
            }


            if (getProperty("checked_snps") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "checked_snps (checkpoints) is not set, we need it to generate a genotype");
                return (null);
            } else {
                this.checkedSNPs = getProperty("checked_snps");
            }

            if (getProperty("check_points") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "check_points (checkpoints) is not set, we need it to generate a report");
                return (null);
            } else {
                this.check_points = getProperty("check_points");
            }

            if (getProperty("existing_matrix") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "existing_matrix not provided, will calclate jaccard indexes for all genotypes");
            } else {
                this.existingMatrix = getProperty("existing_matrix");
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

            if (getProperty("vcftools_version") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "vcftools_version is not set, we need it to call vcftools correctly");
                return (null);
            } else {
                this.vcftoolsVersion = getProperty("vcftools_version");
            }

            if (getProperty("samtools_version") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "samtools_version is not set, we need it to call samtools correctly");
                return (null);
            } else {
                this.samtoolsVersion = getProperty("samtools_version");
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

            if (getProperty("provision_vcfs") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "provision_vcfs flag is not set, will not provision vcf files");
            } else {
                this.provisionVcfs = getProperty("provision_vcfs").toString().isEmpty() ? false : true;
            }

            if (getProperty("manual_output") == null) {
                this.manualOutput = false;
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "No manual output requested, will append a random dir to the output path");
            } else {
                String manualCheck = getProperty("manual_output");
                this.manualOutput = (manualCheck.isEmpty() || manualCheck.equalsIgnoreCase("false")) ? false : true;
            }

            this.vcf_files = new String[this.bam_files.length];

            for (int i = 0; i < this.bam_files.length; i++) {
                //Using first file, try to guess the study name if it was not provided as an argument in .ini file
                if (i == 0 && this.studyName.isEmpty()) {
                    if (bam_files[i].matches("SWID_\\d+_\\D+_\\d+")) {
                        String tempName = bam_files[i].substring(bam_files[i].lastIndexOf("SWID_"));
                        this.studyName = tempName.substring(0, tempName.indexOf("_") - 1);
                    } else {
                        this.studyName = "UNDEF";
                    }
                }

                String basename = this.bam_files[i].substring(this.bam_files[i].lastIndexOf("/") + 1, this.bam_files[i].lastIndexOf(".bam"));
                if (basename.contains(".")) {
                    basename = basename.substring(0, basename.indexOf("."));
                }
                String vcfName = basename + ".snps.raw.vcf";
                this.vcf_files[i] = this.dataDir + vcfName;

                // If we don't have a vcf file we need to set it as an output (Obsolete)
                // vcf files are in the same order as bam files, the decider needs to ensure this

            }
            Log.stdout("Created array of " + this.vcf_files.length + " vcf files of length ");

            //Construct GATK_dirs array with random names for later use with GATK
            this.GATK_dirs = new String[this.bam_files.length];
            int seed = this.makeRandom(8);
            for (int b = 0; b < this.bam_files.length; b++) {
                this.GATK_dirs[b] = this.tmpDir + seed++;
            }
        } catch (Exception e) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, null, e);
        }
        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        try {
            this.dataDir = getProperty("data_dir");
            if (this.dataDir == null || this.dataDir.isEmpty()) {
                this.dataDir = "data/";
            }
            if (!this.dataDir.endsWith("/")) {
                this.dataDir += "/";
            }
            this.addDirectory(this.dataDir);
            this.addDirectory(this.tempDir);
            this.addDirectory(this.dataDir + this.finDir);
            int num_of_files = getProperty("input_files").split(",").length;
            //Make a pool of tmp directories for GATK:
            this.GATK_dirs = new String[num_of_files];
            int seed = this.makeRandom(8);
            for (int b = 0; b < num_of_files; b++) {
                this.GATK_dirs[b] = this.tmpDir + seed++;
                this.addDirectory(GATK_dirs[b]);
            }

        } catch (Exception e) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, null, e);
        }

    }

    @Override
    public void buildWorkflow() {
        try {
            //Need the decider to check the headers of the .bam files for the presence of RG and abscence of empty Fields (Field: )
            Workflow workflow = this.getWorkflow();
            int newVcfs = 0;
            List<Job> prov_jobs = new ArrayList<Job>();
            List<Job> gatk_jobs = new ArrayList<Job>();

            // A small copy job
            Job job_copy = workflow.createBashJob("copy_resources");
            job_copy.setCommand("cp -R " + getWorkflowBaseDir()
                    + "/data/* "
                    + this.dataDir);

            job_copy.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                job_copy.setQueue(this.queue);
            }

            // TODO Entry point for Job batching modifications
            int batchLength = this.bam_files.length / this.batchCount;
            //int batchLength = 15; //DEBUG
            String[] bam_path = new String[this.bam_files.length];

            for (int i = 0; i < this.bam_files.length; i += batchLength) {
                //BATCHING Job1 series
                int stop = i + batchLength >= this.bam_files.length ? this.bam_files.length : i + batchLength;
                prov_jobs.clear();
                
                for (int bj0 = i; bj0 < stop; bj0++) {  
                    
                  Job job_provide = workflow.createBashJob("provide_bams_" + i);
                  job_provide.setCommand("echo Provisioning File in Batch " + i);
                  StringBuilder provCommand = new StringBuilder();

                  SqwFile bamFile = this.createInFile("application/bam", this.bam_files[bj0]);
                  bam_path[bj0] = bamFile.getProvisionedPath();
                  job_provide.addFile(bamFile);

                  provCommand.append("echo Provisioning File ").append(bj0);

                  job_provide.setCommand(provCommand.toString());
                  job_provide.addParent(job_copy);
                  job_provide.setMaxMemory("2000");
                  if (!this.queue.isEmpty()) {
                    job_provide.setQueue(this.queue);
                  }
                  prov_jobs.add(job_provide);
                } //Batching ENDs
                
                Job job_index = workflow.createBashJob("index_bams_" + i);
                StringBuilder idxCommand = new StringBuilder();
                for (int bj1 = i; bj1 < stop; bj1++) {  

                  if (bj1 > i)
                      idxCommand.append(" && ");

                  idxCommand.append(getWorkflowBaseDir())
                            .append("/bin/samtools-")
                            .append(this.samtoolsVersion)
                            .append("/samtools index ")
                            .append(bam_path[bj1]);
               
                } //Batching ENDs
                job_index.setCommand(idxCommand.toString());
                job_index.setMaxMemory("4000");
                if (!this.queue.isEmpty()) {
                    job_index.setQueue(this.queue);
                }
                for (Job pj : prov_jobs) {
                    job_index.addParent(pj);
                }
                
                Job job_gatk = workflow.createBashJob("call_snps_" + i);
                StringBuilder gatkCommand = new StringBuilder();
                //BATCHING Job2 series
                for (int bj2 = i; bj2 < stop; bj2++) {
                    if (bj2 > i) {
                        gatkCommand.append(" && ");
                    }
                    gatkCommand.append(gatk_java)
                            .append(" -Xmx2g -Djava.io.tmpdir=").append(this.GATK_dirs[bj2])
                            .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                            .append("-R ").append(this.genomeFile).append(" ")
                            .append("-T UnifiedGenotyper -I ")
                            .append(bam_path[bj2]).append(" ")
                            .append("-o ").append(this.vcf_files[i]).append(" ")
                            .append("-stand_call_conf ").append(this.stand_call_conf).append(" ")
                            .append("-stand_emit_conf ").append(this.stand_emit_conf).append(" ")
                            .append("-dcov ").append(this.dcov).append(" ")
                            .append("-L ").append(this.checkedSNPs);
                    newVcfs++; // Used only to track the number of newly generated vcf files 
                    if (this.provisionVcfs) {
                        SqwFile vcf_file = this.createOutputFile(this.vcf_files[bj2], "text/vcf-4", this.manualOutput);
                        job_gatk.addFile(vcf_file);
                    }
                } //Batching ENDs
                job_gatk.setCommand(gatkCommand.toString());
                job_gatk.setMaxMemory(getProperty("gatk_memory"));
                if (!this.queue.isEmpty()) {
                    job_gatk.setQueue(this.queue);
                }
                job_gatk.addParent(job_index);
                gatk_jobs.add(job_gatk);

                Job job_gatk2 = workflow.createBashJob("calculate_depth_" + i);
                StringBuilder gatk2Command = new StringBuilder();
                //BATCHING Job3 series
                for (int bj3 = i; bj3 < stop; bj3++) {

                    if (bj3 > i) {
                        gatk2Command.append(" && ");
                    }
                    gatk2Command.append(gatk_java).append(" -Xmx3g -Djava.io.tmpdir=").append(this.GATK_dirs[bj3])
                            .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                            .append("-R ").append(this.genomeFile).append(" ")
                            .append("-T DepthOfCoverage ").append("-I ").append(bam_path[bj3]).append(" ")
                            .append("-o ").append(this.tempDir).append(this.makeBasename(this.bam_files[bj3])).append(" ")
                            .append("-L ").append(this.checkedSNPs);

                }// BATCHING ENDs
                job_gatk2.setCommand(gatk2Command.toString());
                job_gatk2.setMaxMemory(getProperty("gatk_memory"));
                if (!this.queue.isEmpty()) {
                    job_gatk2.setQueue(this.queue);
                }
                job_gatk2.addParent(job_index);

                //Additional step to create depth of coverage data - for individual fingerprint image generation
                Job job_fin = workflow.createBashJob("make_fingerprint_file_" + i);
                StringBuilder finCommand = new StringBuilder();
                //BATCHING Job4 series
                for (int bj4 = i; bj4 < stop; bj4++) {

                    if (bj4 > i) {
                        finCommand.append(" && ");
                    }
                    String basename = this.makeBasename(this.bam_files[bj4]);
                    finCommand.append(getWorkflowBaseDir()).append("/dependencies/create_fin.pl")
                              .append(" --refvcf=").append(this.checkedSNPs)
                              .append(" --genotype=").append(this.vcf_files[bj4])
                              .append(" --coverage=").append(this.tempDir).append(basename).append(".sample_interval_summary")
                              .append(" --datadir=").append(this.tempDir)
                              .append(" --outdir=").append(this.dataDir).append(this.finDir)
                              .append(" --basename=").append(basename);
                } //BATCHING ENDs
                job_fin.setCommand(finCommand.toString());
                job_fin.setMaxMemory("4000");
                job_fin.addParent(job_gatk2);
                gatk_jobs.add(job_fin);
            } // END OF job Batching loop

            // We don't need to continue if there are no new vcf files to generate
            if (newVcfs == 0) {
                Logger.getLogger(getClass().getName()).log(Level.INFO, "There are no new genotypes to generate, to avoid duplicate calculation the workflow will terminate");
                System.exit(0);
            }

            // At his point we need to stage, bgzip and tabix all vcf files we would need to create a jaccard matrix          
            Job job_vcfprep = workflow.createBashJob("prepare_vcfs");
            job_vcfprep.setCommand(getWorkflowBaseDir() + "/dependencies/prepare_vcfs.pl "
                    + "--datadir=" + this.dataDir + " "
                    + "--tabix=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix "
                    + "--bgzip=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip");
            job_vcfprep.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                job_vcfprep.setQueue(this.queue);
            }

            //SET PARENTS FOR VCFPREP JOB
            for (Job parent : gatk_jobs) {
                job_vcfprep.addParent(parent);
            }

            List<Job> upstream_jobs = new ArrayList<Job>();
            upstream_jobs.add(job_vcfprep);


            // We use simple reasoning here - if existingMatrix is unavailable, we will need to run a loooong job
            // of calculting jaccard indexes if the list of vcf files is long. So, we will break the list into smaller chunks
            // And use the results as exisitngMatrix in the next step Likely to run only once
            // (when running this workflow on a study for the very first time)

            if (this.existingMatrix.isEmpty() && this.vcf_files.length > this.jChunkSize) {
                int chunkMultiplier = 1;

                StringBuilder chunkedResults = new StringBuilder();
                List<Integer> vcf_chunks = new ArrayList<Integer>(this.vcf_files.length);

                for (int vc = 1; vc < this.vcf_files.length; vc++) {
                    if (vc < this.jChunkSize * chunkMultiplier) {
                        continue;
                    }
                    vcf_chunks.add(new Integer(1));
                    chunkMultiplier += 1;
                }
                vcf_chunks.add(new Integer(1));
                int resultID = 1;

                // Combine chunks and run the jobs, registering results for later use
                for (int c = 0; c < vcf_chunks.size(); c++) {
                    for (int cc = c; cc < vcf_chunks.size(); cc++) {
                        if (c == cc) {
                            continue;
                        }

                        Job job_list_writer = workflow.createBashJob("make_list" + resultID);
                        String chunkList = "jaccard.chunk." + resultID + ".list";
                        job_list_writer.setCommand(getWorkflowBaseDir() + "/dependencies/write_list.pl "
                                + "--datadir=" + this.dataDir + " "
                                + "--segments=\"" + c * this.jChunkSize + ":" + ((c + 1) * this.jChunkSize - 1) + ","
                                + cc * this.jChunkSize + ":" + (cc + 1) * this.jChunkSize + "\" "
                                + "> " + this.dataDir + chunkList);

                        job_list_writer.setMaxMemory("2000");
                        if (!this.queue.isEmpty()) {
                            job_list_writer.setQueue(this.queue);
                        }

                        job_list_writer.addParent(job_vcfprep);

                        Job job_jchunk = workflow.createBashJob("make_matrix_" + resultID);
                        String chunkName = "jaccard.chunk." + resultID + ".csv";
                        String sep = chunkedResults.toString().isEmpty() ? "" : ",";
                        chunkedResults.append(sep).append(this.dataDir).append(chunkName);

                        resultID += 1;

                        job_jchunk.setCommand(getWorkflowBaseDir() + "/dependencies/jaccard_coeff.matrix.pl "
                                + "--list=" + this.dataDir + chunkList + " "
                                + "--vcf_compare=" + getWorkflowBaseDir() + "/bin/vcftools_" + this.vcftoolsVersion + "/bin/vcf-compare "
                                + "--datadir=" + this.dataDir + " "
                                + "--tabix=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + " "
                                + "--studyname=" + this.studyName + " "
                                + "> " + this.dataDir + chunkName);

                        job_jchunk.setMaxMemory("2000");
                        if (!this.queue.isEmpty()) {
                            job_jchunk.setQueue(this.queue);
                        }

                        job_jchunk.addParent(job_list_writer);
                        upstream_jobs.add(job_jchunk);
                    }
                }
                // The result of the last job will be existingMatrix
                this.existingMatrix = chunkedResults.toString();
            }


            // Next job is jaccard_coeff.matrix script, need to create|update giant matrix of jaccard indexes
            // using vcftools, colors should not be assigned at his step. Will be final jaccard index job          
            Job job_list_writer2 = workflow.createBashJob("make_Final_list");
            String chunkList = "jaccard.chunks.list";
            job_list_writer2.setCommand(getWorkflowBaseDir() + "/dependencies/write_list.pl "
                    + "--datadir=" + this.dataDir + " "
                    + "> " + this.dataDir + chunkList);

            job_list_writer2.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                job_list_writer2.setQueue(this.queue);
            }
            job_list_writer2.addParent(job_vcfprep);

            SqwFile matrix = this.createOutputFile(this.dataDir + this.studyName + "_jaccard.matrix.csv", "text/plain", this.manualOutput);
            Job job_jaccard = workflow.createBashJob("make_matrix");
            job_jaccard.setCommand(getWorkflowBaseDir() + "/dependencies/jaccard_coeff.matrix.pl "
                    + "--list=" + this.dataDir + chunkList + " "
                    + "--vcf_compare=" + getWorkflowBaseDir() + "/bin/vcftools_" + this.vcftoolsVersion + "/bin/vcf-compare "
                    + "--datadir=" + this.dataDir + " "
                    + "--tabix=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + " "
                    + "--studyname=" + this.studyName + " "
                    + "> " + matrix.getSourcePath());
            if (!this.existingMatrix.isEmpty()) {
                job_jaccard.getCommand().addArgument("--existing_matrix=" + this.existingMatrix);
            }

            job_jaccard.setMaxMemory("3000");
            if (!this.queue.isEmpty()) {
                job_jaccard.setQueue(this.queue);
            }
            job_jaccard.addFile(matrix);
            job_jaccard.addParent(job_list_writer2);
            for (Job parent : upstream_jobs) {
                job_jaccard.addParent(parent);
            }

            // Images generated here: matrix slicing/color assignment, flagging suspicious files, wrapper for R script
            Job make_pics = workflow.createBashJob("make_report");
            make_pics.setCommand("perl " + getWorkflowBaseDir() + "/dependencies/make_report.pl "
                    + "--matrix=" + matrix.getSourcePath() + " "
                    + "--refsnps=" + this.check_points + " "
                    + "--tempdir=" + this.dataDir + this.finDir + " "
                    + "--datadir=" + this.dataDir + " "
                    + "--studyname=" + this.studyName + " "
                    + "> " + this.dataDir + "index.html");
            make_pics.addParent(job_copy);
            make_pics.addParent(job_jaccard);
            make_pics.setMaxMemory("4000");
            if (!this.queue.isEmpty()) {
                make_pics.setQueue(this.queue);
            }

            // Zip finfiles and similarity matrix for customization in the webtool
            Job zip_fins = workflow.createBashJob("zip_finfiles");
            zip_fins.setCommand("zip -r " + this.dataDir + "customize.me.zip "
                    + this.dataDir + "finfiles "
                    + matrix.getSourcePath());
            zip_fins.addParent(make_pics);
            zip_fins.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                zip_fins.setQueue(this.queue);
            }

            // Zip everything into a report bundle for provisioning
            String report_name = this.reportName + "." + this.studyName;
            SqwFile report_file = this.createOutputFile(this.dataDir + report_name + ".report.zip", "application/zip-report-bundle", manualOutput);
            Job zip_report = workflow.createBashJob("zip_everything");
            zip_report.setCommand("zip -r " + report_file.getSourcePath() + " "
                    + this.dataDir + "*.png "
                    + this.dataDir + "images "
                    + this.dataDir + "*_genotype_report*.csv "
                    + this.dataDir + "*_similarity_matrix*.csv "
                    + matrix.getSourcePath() + " "
                    + this.dataDir + "customize.me.zip "
                    + this.dataDir + "*.html");
            zip_report.addParent(zip_fins);
            zip_report.setMaxMemory("2000");

            zip_report.addFile(report_file);

            if (!this.queue.isEmpty()) {
                zip_report.setQueue(this.queue);
            }

            // Next job is report-emailing script (check html for flagged samples, send email(s) to interested people)
            if (!this.watchersList.isEmpty()) {
                Job job_alert = workflow.createBashJob("send_alert");
                job_alert.setCommand("perl " + getWorkflowBaseDir() + "/dependencies/plotReporter.pl "
                        + "--watchers=" + this.watchersList + " "
                        + "--report=" + this.dataDir + "index.html "
                        + "--studyname=" + this.studyName + " "
                        + "--bundle=" + report_file.getOutputPath());
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

    private SqwFile createInFile(String meta, String source) {
        SqwFile file = new SqwFile();
        file.setType(meta);
        file.setSourcePath(source);
        file.setIsInput(true);
        return file;
    }

    private int makeRandom(int digits) {
        Random rnd = new Random();
        String randomString = "1"; // start each string from 1 to prevent zero at the beginning 
        for (int i = 1; i < digits; i++) {
            randomString += rnd.nextInt(9);
        }
        return Integer.parseInt(randomString);
    }

    private String makeBasename(String name) {
        String basename = name.substring(name.lastIndexOf("/") + 1, name.lastIndexOf(".bam"));
        if (basename.contains(".")) {
            basename = basename.substring(0, basename.indexOf("."));
        }
        return basename;
    }
}
