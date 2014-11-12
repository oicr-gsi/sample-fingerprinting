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
 * Workflow. For pairing with SampleFingerprinting 2.0 A lighter (split) version
 * of SampleFingerprinting 1.1
 * 
 * @author pruzanov
 */
public class FingerprintCollectorWorkflow extends OicrWorkflow {

    private String[] bamFiles;
    private String[] vcfFiles;
    private String[] gatkDirs;
    private String studyName = "";
    //private String watchersList = "";
    private String dataDir;
    private String tempDir = "tempfiles/";
    private String gatkPrefix = "./";
    //Additional one for GATK:
    private String gatkTmp = "temp";
    private String finDir = "finfiles/";
    private String gatkVersion;
    private String tabixVersion;
    private String vcftoolsVersion;
    private String samtoolsVersion;
    private String genomeFile;
    private String checkedSNPs;
    private String queue;
    private String reportName = "sample_fingerprint";
    private String gatkJava;
    //GATK parameters
    private String standCallConf = "50.0";
    private String standEmitConf = "10.0";
    private String dcov = "200";
    //private final int jChunkSize = 50; // Maximum allowed number of vcf files when jaccard_indexing step doesn't fork into multiple sub-jobs
    private final int batchCount = 100; // Use for job batching, this many jobs
    private boolean manualOutput;
    private boolean provisionVcfs;

    //Misc
    private static final int RANDOM_SEED = 8;
    private static final int RANDOM_PICK = 9;

    @Override
    public Map<String, SqwFile> setupFiles() {
        // Set up reference, bam and vcf files here
        try {
           if (getProperty("input_files") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "input_files is not set, we need at least one bam file");
                return (null);
            } else {
                this.bamFiles = getProperty("input_files").split(",");
            }

            if (getProperty("gatk_java") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "gatk_java is not set, we need it to run GATK since it requires 1.7 java and we cannot rely on the defaul");
                return (null);
            } else {
                this.gatkJava = getProperty("gatk_java");
            }

            if (getProperty("genome_file") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "genome_file (reference assembly fasta) is not set, we need it to generate a genotype");
                return (null);
            } else {
                this.genomeFile = getProperty("genome_file");
            }


            if (getProperty("checked_snps") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "checked_snps (checkpoints) is not set, we need it to generate a genotype");
                return (null);
            } else {
                this.checkedSNPs = getProperty("checked_snps");
            }

            /*
            if (getProperty("existing_matrix") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "existing_matrix not provided, will calclate jaccard indexes for all genotypes");
            } else {
                this.existingMatrix = getProperty("existing_matrix");
            }

            if (getProperty("watchers_list") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "watchers list not provided, will not send alerts");
            } else {
                this.watchersList = getProperty("watchers_list");
            }
            */

            this.standCallConf = getOptionalProperty("stand_call_conf", "50.0");
            this.standEmitConf = getOptionalProperty("stand_emit_conf", "10.0");
            this.dcov            = getOptionalProperty("dcov", "200");

            if (getProperty("gatk_version") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "gatk_version is not set, we need it to call GATK correctly");
                return (null);
            } else {
                this.gatkVersion = getProperty("gatk_version");
            }

            if (getProperty("tabix_version") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "tabix_version is not set, we need it to call tabix correctly");
                return (null);
            } else {
                this.tabixVersion = getProperty("tabix_version");
            }

            if (getProperty("vcftools_version") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "vcftools_version is not set, we need it to call vcftools correctly");
                return (null);
            } else {
                this.vcftoolsVersion = getProperty("vcftools_version");
            }

            if (getProperty("samtools_version") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "samtools_version is not set, we need it to call samtools correctly");
                return (null);
            } else {
                this.samtoolsVersion = getProperty("samtools_version");
            }

            if (getProperty("queue") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "Queue not set, most likely will run as default queue");
                this.queue = "";
            } else {
                this.queue = getProperty("queue");
            }

            if (getProperty("study_name") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "Study name isn't not set, will try to extract it from file names");
                this.studyName = "";
            } else {
                this.studyName = getProperty("study_name");
            }

            if (getProperty("provision_vcfs") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "provision_vcfs flag is not set, will not provision vcf files");
            } else {
                this.provisionVcfs = getProperty("provision_vcfs").toString().isEmpty() ? false : true;
            }

            if (getProperty("manual_output") == null) {
                this.manualOutput = false;
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "No manual output requested, will append a random dir to the output path");
            } else {
                String manualCheck = getProperty("manual_output");
                this.manualOutput = (manualCheck.isEmpty() || manualCheck.equalsIgnoreCase("false")) ? false : true;
            }

            this.vcfFiles = new String[this.bamFiles.length];

            for (int i = 0; i < this.bamFiles.length; i++) {
                //Using first file, try to guess the study name if it was not provided as an argument in .ini file
                if (i == 0 && this.studyName.isEmpty()) {
                    if (bamFiles[i].matches("SWID_\\d+_\\D+_\\d+")) {
                        String tempName = bamFiles[i].substring(bamFiles[i].lastIndexOf("SWID_"));
                        this.studyName = tempName.substring(0, tempName.indexOf("_") - 1);
                    } else {
                        this.studyName = "UNDEF";
                    }
                }

                String basename = this.bamFiles[i].substring(this.bamFiles[i].lastIndexOf("/") + 1, this.bamFiles[i].lastIndexOf(".bam"));
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
            Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, null, e);
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
            //Setup gatk prefix
            this.gatkPrefix = getOptionalProperty("gatk_prefix", "./");
            if (!this.gatkPrefix.endsWith("/")) {
                     this.gatkPrefix += "/";
            }

            this.addDirectory(this.dataDir);
            this.addDirectory(this.tempDir);
            this.addDirectory(this.dataDir + this.finDir);
            int numberOfFiles = getProperty("input_files").split(",").length;
            //Make a pool of tmp directories for GATK:
            this.gatkDirs = new String[numberOfFiles];
            int seed = this.makeRandom(RANDOM_SEED);
            for (int b = 0; b < numberOfFiles; b++) {
                this.gatkDirs[b] = this.gatkPrefix + this.gatkTmp + seed++;
                this.addDirectory(gatkDirs[b]);
            }

        } catch (Exception e) {
            Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, null, e);
        }

    }

    @Override
    public void buildWorkflow() {
        try {
            //Need the decider to check the headers of the .bam files for the presence of RG and abscence of empty Fields (Field: )
            Workflow workflow = this.getWorkflow();
            int newVcfs = 0;
            List<Job> provisionJobs = new ArrayList<Job>();
            List<Job> gatkJobs = new ArrayList<Job>();

            // TODO Entry point for Job batching modifications
            int batchLength = this.bamFiles.length / this.batchCount;
            String[] bamfilePaths = new String[this.bamFiles.length];

            for (int i = 0; i < this.bamFiles.length; i += batchLength) {
                //BATCHING Job1 series
                int stop = i + batchLength >= this.bamFiles.length ? this.bamFiles.length : i + batchLength;
                provisionJobs.clear();

                for (int bj0 = i; bj0 < stop; bj0++) {

                  Job jobProvide = workflow.createBashJob("provide_bams_" + bj0);
                  jobProvide.setCommand("echo Provisioning File in Batch " + i);
                  StringBuilder provCommand = new StringBuilder();

                  SqwFile bamFile = this.createInFile("application/bam", this.bamFiles[bj0]);
                  bamfilePaths[bj0] = bamFile.getProvisionedPath();
                  jobProvide.addFile(bamFile);

                  provCommand.append("echo Provisioning File ").append(bj0);

                  jobProvide.setCommand(provCommand.toString());
                  jobProvide.setMaxMemory("2000");
                  if (!this.queue.isEmpty()) {
                    jobProvide.setQueue(this.queue);
                  }
                  provisionJobs.add(jobProvide);
                } //Batching ENDs

                Job jobIndex = workflow.createBashJob("index_bams_" + i);
                StringBuilder idxCommand = new StringBuilder();
                for (int bj1 = i; bj1 < stop; bj1++) {

                  if (bj1 > i) {
                      idxCommand.append(" && ");
                  }
                  idxCommand.append(getWorkflowBaseDir())
                            .append("/bin/samtools-")
                            .append(this.samtoolsVersion)
                            .append("/samtools index ")
                            .append(bamfilePaths[bj1]);

                } //Batching ENDs
                jobIndex.setCommand(idxCommand.toString());
                jobIndex.setMaxMemory("4000");
                if (!this.queue.isEmpty()) {
                    jobIndex.setQueue(this.queue);
                }
                for (Job pj : provisionJobs) {
                    jobIndex.addParent(pj);
                }

                Job jobGATK = workflow.createBashJob("call_snps_" + i);
                StringBuilder gatkCommand = new StringBuilder();
                //BATCHING Job2 series
                for (int bj2 = i; bj2 < stop; bj2++) {
                    if (bj2 > i) {
                        gatkCommand.append(" && ");
                    }
                    gatkCommand.append(gatkJava)
                            .append(" -Xmx2g -Djava.io.tmpdir=").append(this.gatkDirs[bj2])
                            .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                            .append("-R ").append(this.genomeFile).append(" ")
                            .append("-T UnifiedGenotyper -I ")
                            .append(bamfilePaths[bj2]).append(" ")
                            .append("-o ").append(this.vcfFiles[bj2]).append(" ")
                            .append("-stand_call_conf ").append(this.standCallConf).append(" ")
                            .append("-stand_emit_conf ").append(this.standEmitConf).append(" ")
                            .append("-dcov ").append(this.dcov).append(" ")
                            .append("-L ").append(this.checkedSNPs);
                    newVcfs++; // Used only to track the number of newly generated vcf files
                    if (this.provisionVcfs) {
                        SqwFile vcfFile = this.createOutputFile(this.vcfFiles[bj2], "text/vcf-4", this.manualOutput);
                        jobGATK.addFile(vcfFile);
                    }
                } //Batching ENDs
                jobGATK.setCommand(gatkCommand.toString());
                jobGATK.setMaxMemory(getProperty("gatk_memory"));
                if (!this.queue.isEmpty()) {
                    jobGATK.setQueue(this.queue);
                }
                jobGATK.addParent(jobIndex);
                gatkJobs.add(jobGATK);

                Job jobGATK2 = workflow.createBashJob("calculate_depth_" + i);
                StringBuilder gatk2Command = new StringBuilder();
                //BATCHING Job3 series
                for (int bj3 = i; bj3 < stop; bj3++) {

                    if (bj3 > i) {
                        gatk2Command.append(" && ");
                    }
                    gatk2Command.append(gatkJava).append(" -Xmx3g -Djava.io.tmpdir=").append(this.gatkDirs[bj3])
                            .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                            .append("-R ").append(this.genomeFile).append(" ")
                            .append("-T DepthOfCoverage ").append("-I ").append(bamfilePaths[bj3]).append(" ")
                            .append("-o ").append(this.tempDir).append(this.makeBasename(this.bamFiles[bj3])).append(" ")
                            .append("-L ").append(this.checkedSNPs);
                // BATCHING ENDs
                }
                jobGATK2.setCommand(gatk2Command.toString());
                jobGATK2.setMaxMemory(getProperty("gatk_memory"));
                if (!this.queue.isEmpty()) {
                    jobGATK2.setQueue(this.queue);
                }
                jobGATK2.addParent(jobIndex);
                //Additional step to create depth of coverage data - for individual fingerprint image generation
                Job jobFin = workflow.createBashJob("make_fingerprint_file_" + i);
                StringBuilder finCommand = new StringBuilder();
                //BATCHING Job4 series
                for (int bj4 = i; bj4 < stop; bj4++) {

                    if (bj4 > i) {
                        finCommand.append(" && ");
                    }
                    String basename = this.makeBasename(this.bamFiles[bj4]);
                    finCommand.append(getWorkflowBaseDir()).append("/dependencies/create_fin.pl")
                              .append(" --refvcf=").append(this.checkedSNPs)
                              .append(" --genotype=").append(this.vcfFiles[bj4])
                              .append(" --coverage=").append(this.tempDir).append(basename).append(".sample_interval_summary")
                              .append(" --datadir=").append(this.tempDir)
                              .append(" --outdir=").append(this.dataDir).append(this.finDir)
                              .append(" --basename=").append(basename);
                } //BATCHING ENDs
                jobFin.setCommand(finCommand.toString());
                jobFin.setMaxMemory("4000");
                jobFin.addParent(jobGATK2);
                gatkJobs.add(jobFin);
            } // END OF job Batching loop

            // We don't need to continue if there are no new vcf files to generate
            if (newVcfs == 0) {
                Logger.getLogger(getClass().getName()).log(Level.INFO, "There are no new genotypes to generate, to avoid duplicate calculation the workflow will terminate");
                throw new RuntimeException("There is not enough of new data, terminating...");
            }

            // At his point we need to stage, bgzip and tabix all vcf files we would need to create a jaccard matrix
            Job jobVcfPrep = workflow.createBashJob("prepare_vcfs");
            jobVcfPrep.setCommand(getWorkflowBaseDir() + "/dependencies/prepare_vcfs.pl "
                    + "--datadir=" + this.dataDir + " "
                    + "--tabix=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/tabix "
                    + "--bgzip=" + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion + "/bgzip");
            jobVcfPrep.setMaxMemory("2000");
            if (!this.queue.isEmpty()) {
                jobVcfPrep.setQueue(this.queue);
            }

            //SET PARENTS FOR VCFPREP JOB
            for (Job parent : gatkJobs) {
                jobVcfPrep.addParent(parent);
            }

            // TODO Need to set up a provsioning job that may be preceeded by
            // zipping all needed files into a single archive

            /*
             * ======== Next, SampleFingerprinting 2.0 should kick in ==========
             */


        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }

    }

    private SqwFile createInFile(final String meta, final String source) {
        SqwFile file = new SqwFile();
        file.setType(meta);
        file.setSourcePath(source);
        file.setIsInput(true);
        return file;
    }

    private int makeRandom(final int digits) {
        Random rnd = new Random();
        StringBuilder sb = new StringBuilder();
        sb.append("1"); // start each string from 1 to prevent zero at the beginning
        for (int i = 1; i < digits; i++) {
            sb.append(rnd.nextInt(RANDOM_PICK));
        }
        return Integer.parseInt(sb.toString());
    }

    private String makeBasename(final String name) {
        String basename = name.substring(name.lastIndexOf("/") + 1, name.lastIndexOf(".bam"));
        if (basename.contains(".")) {
            basename = basename.substring(0, basename.indexOf("."));
        }
        return basename;
    }
}
