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
 * SampleFingerprinting v 2.0 - workflow split from original (1.1) version and
 * working with Fingerprint Collector workflow
 * 
 * @author pruzanov
 */
public class SampleFingerprintingWorkflow extends OicrWorkflow {

  //  private String[] bam_files;
  //  private String[] genotypes;
    private String[] vcfFiles;
    private String[] finFiles;
   // private String[] GATK_dirs;
    private String existingMatrix = "";
    private String studyName = "";
    private String watchersList = "";
    private String dataDir;
    //private String tempDir = "tempfiles/";
 //   private String gatkPrefix = "./";
    //Additional one for GATK:
    //private String gatkTmp = "temp";
    private final String finDir = "finfiles/";
   // private String gatkVersion;
    private String tabixVersion;
    private String vcftoolsVersion;
   // private String samtoolsVersion;
   // private String genomeFile;
   // private String checkedSNPs;
    private String check_points;
    private String queue;
    private final String reportName = "sample_fingerprint";
   // private String gatk_java;
    //GATK parameters
  //  private String stand_call_conf = "50.0";
  //  private String stand_emit_conf = "10.0";
  //  private String dcov = "200";
    private final int jChunkSize = 50; // Maximum allowed number of vcf files when jaccard_indexing step doesn't fork into multiple sub-jobs
  //  private final int batchCount = 100; // Use for job batching, this many jobs 
    private boolean manualOutput;
    
    private final static String VCF_EXT   = ".snps.raw.vcf.gz";
    private final static String TABIX_EXT = ".snps.raw.vcf.gz.tbi";
    private final static String FIN_EXT   = ".fin";

    @Override
    public Map<String, SqwFile> setupFiles() {
        // Set up reference, bam and vcf files here     
        try {
           if (getProperty("input_files") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "input_files is not set, we need at least one bam file");
                return (null);
            } else {
                this.vcfFiles = getProperty("input_files").split(",");
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

            if (getProperty("manual_output") == null) {
                this.manualOutput = false;
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "No manual output requested, will append a random dir to the output path");
            } else {
                String manualCheck = getOptionalProperty("manual_output", "false");
                try {
                    this.manualOutput = Boolean.valueOf(manualCheck);
                } catch (NumberFormatException e) {
                    this.manualOutput = false;
                }
            }

            // Set up all inputs, assume that all vcf.gz, vcf.gz.tbi and .fin files sit in triplets in the same directories
            for (int i = 0; i < this.vcfFiles.length; i++) {
                //Using first file, try to guess the study name if it was not provided as an argument in .ini file
                if (i == 0 && this.studyName.isEmpty()) {
                    if (vcfFiles[i].matches("SWID_\\d+_\\D+_\\d+")) {
                        String tempName = vcfFiles[i].substring(vcfFiles[i].lastIndexOf("SWID_"));
                        this.studyName = tempName.substring(0, tempName.indexOf("_") - 1);
                    } else {
                        this.studyName = "UNDEF";
                    }
                }

                String basename = this.vcfFiles[i].substring(this.vcfFiles[i].lastIndexOf("/") + 1, this.vcfFiles[i].lastIndexOf(".vcf.gz"));
                if (basename.contains(".")) {
                    basename = basename.substring(0, basename.indexOf("."));
                }
                String vcfName = basename + ".snps.raw.vcf";
                this.vcfFiles[i] = this.dataDir + vcfName;

                // If we don't have a vcf file we need to set it as an output (Obsolete)
                // vcf files are in the same order as bam files, the decider needs to ensure this

            }
            Log.stdout("Created array of " + this.vcfFiles.length + " vcf files of length ");
        } catch (Exception e) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, null, e);
        }
        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        try {
            //Setup data dir
            this.dataDir = getOptionalProperty("data_dir", "data/");
            if (!this.dataDir.endsWith("/")) {
                this.dataDir += "/";
            }
            this.addDirectory(this.dataDir);
            this.addDirectory(this.finDir);
        } catch (Exception e) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, null, e);
        }

    }

    @Override
    public void buildWorkflow() {
        try {
            //Need the decider to check the headers of the .bam files for the presence of RG and abscence of empty Fields (Field: )
            Workflow workflow = this.getWorkflow();
           // int newVcfs = 0;
            List<Job> prov_jobs = new ArrayList<Job>();
            List<Job> gatk_jobs = new ArrayList<Job>();


            if (this.existingMatrix.isEmpty() && this.vcfFiles.length > this.jChunkSize) {

                StringBuilder chunkedResults = new StringBuilder();
                int vcf_chunks = (int) Math.ceil((double) this.vcfFiles.length / (double) this.jChunkSize);
                int resultID = 1;

                // Combine chunks and run the jobs, registering results for later use
                for (int c = 0; c < vcf_chunks; c++) {
                    for (int cc = c; cc < vcf_chunks; cc++) {

                        Job job_list_writer = workflow.createBashJob("make_list" + resultID);
                        String chunkList = "jaccard.chunk." + resultID + ".list";
                        job_list_writer.setCommand(getWorkflowBaseDir() + "/dependencies/write_list.pl "
                                + "--datadir=" + this.dataDir + " "
                                + "--segments=\"" + c * this.jChunkSize + ":" + ((c + 1) * this.jChunkSize - 1) + ","
                                + cc * this.jChunkSize + ":" + ((cc + 1) * this.jChunkSize - 1) + "\" "
                                + "> " + this.dataDir + chunkList);

                        job_list_writer.setMaxMemory("2000");
                        if (!this.queue.isEmpty()) {
                            job_list_writer.setQueue(this.queue);
                        }

                       // job_list_writer.addParent(job_vcfprep);

                        Job job_jchunk = workflow.createBashJob("make_matrix_" + resultID);
                        String chunkName = "jaccard.chunk." + resultID + ".csv";
                        String sep = chunkedResults.toString().isEmpty() ? "" : ",";
                        chunkedResults.append(sep).append(this.dataDir).append(chunkName);

                        resultID++;

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
                        //upstream_jobs.add(job_jchunk);
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
                    + "--datadir=" + this.dataDir
                    + " > " + this.dataDir + chunkList);

            job_list_writer2.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                job_list_writer2.setQueue(this.queue);
            }
            //job_list_writer2.addParent(job_vcfprep);

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

            // Images generated here: matrix slicing/color assignment, flagging suspicious files, wrapper for R script
            Job make_pics = workflow.createBashJob("make_report");
            make_pics.setCommand("perl " + getWorkflowBaseDir() + "/dependencies/make_report.pl "
                    + "--matrix=" + matrix.getSourcePath() + " "
                    + "--refsnps=" + this.check_points + " "
                    + "--tempdir=" + this.dataDir + this.finDir + " "
                    + "--datadir=" + this.dataDir + " "
                    + "--studyname=" + this.studyName + " "
                    + "> " + this.dataDir + "index.html");
            //make_pics.addParent(job_copy);
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
        StringBuilder sb = new StringBuilder();
        sb.append("1"); // start each string from 1 to prevent zero at the beginning 
        for (int i = 1; i < digits; i++) {
            sb.append(rnd.nextInt(9));
        }
        return Integer.parseInt(sb.toString());
    }

    private String makeBasename(String name) {
        String basename = name.substring(name.lastIndexOf("/") + 1, name.lastIndexOf(".bam"));
        if (basename.contains(".")) {
            basename = basename.substring(0, basename.indexOf("."));
        }
        return basename;
    }
}
