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
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;

/*
 A workaround tested with Bowtie2Alignemr workflow that allows provisioning multiple files into a single dir:
 //DEBUG
 file.setOutputPath("provisionfiles/CUSTOMDIR/" + this.input_fastq1[fileIndex].substring(
 this.input_fastq1[fileIndex].lastIndexOf("/") + 1, this.input_fastq1[fileIndex].length()));
 file.getOutputPath() //GET the path in workflow
 */
/**
 * @Description For pairing with SampleFingerprinting 2.0 A lighter (split)
 * version of SampleFingerprinting 1.1
 *
 * @author pruzanov
 */
public class FingerprintCollectorWorkflow extends OicrWorkflow {

    private String[] bamFiles;
    private String[] vcfFiles;
    private String[] gatkDirs;
    private String studyName = "";
    private String dataDir;
    private final String tempDir = "tempfiles/";
    private final String gatkTmp = "temp";
    private final String finDir = "finfiles/";
    //private final int batchCount = 100; // Use for job batching, this many jobs
    private String gatkPrefix = "./";

    //Additional one for GATK:
    private String gatkVersion;
    private String tabixVersion;
    private String samtoolsVersion;
    private String genomeFile;
    private String checkedSNPs;
    private String queue;
    private String gatkJava;
    private String picardVersion;

    //GATK parameters
    private String standCallConf;
    private String standEmitConf;
    private String dcov;
    private boolean manualOutput;

    //Misc
    private static final int RANDOM_SEED = 8;
    private static final int RANDOM_PICK = 9;
    private static final String STCALLCONF_DEFAULT = "50.0";
    private static final String STEMIT_DEFAULT = "10.0";
    private static final String DCOV_DEFAULT = "200";
    private static final String SNPFILE_SUFFIX = ".snps.raw.vcf";
    private List<String> baseNames;
    private String bundledJRE;
    private String[] rgDetails;
    private boolean doBamFix = false;

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

            this.standCallConf = getOptionalProperty("stand_call_conf", STCALLCONF_DEFAULT);
            this.standEmitConf = getOptionalProperty("stand_emit_conf", STEMIT_DEFAULT);
            this.dcov = getOptionalProperty("dcov", DCOV_DEFAULT);
            this.bundledJRE = getProperty("bundled_jre");
            this.rgDetails = getProperty("rg_details").split(",");
            this.doBamFix = Boolean.valueOf(getOptionalProperty("preprocess_bam", "false"));

            if (this.bamFiles.length != this.rgDetails.length) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "Don't have the same number of Misc data items for bam files");
                return (null);
            }

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

            if (getProperty("picard_version") == null) {
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.SEVERE, "picard_version is not set, we need it to call picard correctly");
                return (null);
            } else {
                this.picardVersion = getProperty("picard_version");
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

            if (getProperty("manual_output") == null) {
                this.manualOutput = false;
                Logger.getLogger(FingerprintCollectorWorkflow.class.getName()).log(Level.WARNING, "No manual output requested, will append a random dir to the output path");
            } else {
                String manualCheck = getProperty("manual_output");
                this.manualOutput = (manualCheck.isEmpty() || manualCheck.equalsIgnoreCase("false")) ? false : true;
            }

            this.vcfFiles = new String[this.bamFiles.length];
            this.baseNames = new ArrayList<String>();
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

                String vcfName = basename + SNPFILE_SUFFIX;
                this.baseNames.add(basename);
                this.vcfFiles[i] = this.dataDir + vcfName;

                // If we don't have a vcf file we need to set it as an output (Obsolete)
                // vcf files are in the same order as bam files, the decider needs to ensure this
            }
            Log.stdout("Created array of " + this.vcfFiles.length + " vcf files");
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

            List<Job> gatkJobs = new ArrayList<Job>();

            String[] orderedBamfilePaths = new String[this.bamFiles.length];
            String[] bamfilePathsRG = new String[this.bamFiles.length];
            String[] gatkInputs = new String[this.bamFiles.length];

            for (int i = 0; i < this.bamFiles.length; i++) {

                SqwFile bamFile = this.createInFile("application/bam", this.bamFiles[i]);
                bamFile.setIsInput(true);
                List<Job> fixerJobs = new ArrayList<Job>();

                if (this.doBamFix) {
                    Log.stdout("Bam pre-processing requested, will re-order and add Read Groups");
                    orderedBamfilePaths[i] = this.dataDir + this.makeBasename(bamFile.getProvisionedPath()) + ".ordered.bam";
                    bamfilePathsRG[i] = this.dataDir + this.makeBasename(bamFile.getProvisionedPath()) + ".rg.bam";
                    gatkInputs[i] = bamfilePathsRG[i];

                    Job jobReorder = workflow.createBashJob("reorder_bams_" + i);
                    jobReorder.addFile(bamFile);
                    jobReorder.setCommand(getWorkflowBaseDir() + "/bin/" + this.bundledJRE + "/bin/java"
                            + " -Xmx4000M"
                            + " -jar " + getWorkflowBaseDir() + "/bin/picard-tools-" + this.picardVersion + "/ReorderSam.jar"
                            + " INPUT=" + bamFile.getProvisionedPath()
                            + " OUTPUT=" + orderedBamfilePaths[i]
                            + " REFERENCE=" + this.genomeFile);

                    jobReorder.setMaxMemory("7000");

                    if (!this.queue.isEmpty()) {
                        jobReorder.setQueue(this.queue);
                    }

                    Job jobIndex = workflow.createBashJob("index_bams_" + i);

                    String[] rgItems = this.rgDetails[i].split(":");
                    if (rgItems.length != 5) {
                        throw new RuntimeException("RG data entry corrupted for " + this.makeBasename(orderedBamfilePaths[i]) + ".bam, terminating...");
                    }
                    jobIndex.setCommand(getWorkflowBaseDir() + "/bin/" + this.bundledJRE + "/bin/java"
                            + " -Xmx" + this.getProperty("picard_memory") + "M"
                            + " -jar " + getWorkflowBaseDir() + "/bin/picard-tools-" + this.picardVersion + "/AddOrReplaceReadGroups.jar"
                            + " INPUT=" + orderedBamfilePaths[i]
                            + " OUTPUT=" + bamfilePathsRG[i]
                            + " RGID=" + rgItems[0]
                            + " RGLB=" + rgItems[1]
                            + " RGPL=" + rgItems[2]
                            + " RGSM=" + rgItems[3]
                            + " RGPU=" + rgItems[4]
                            + " CREATE_INDEX=TRUE"
                            + " SORT_ORDER=coordinate");

                    jobIndex.setMaxMemory("10000");
                    jobIndex.addParent(jobReorder);
                    if (!this.queue.isEmpty()) {
                        jobIndex.setQueue(this.queue);
                    }
                    fixerJobs.add(jobIndex);
                    
                } else {
                    Log.stdout("Assuming that bam files are ordered and have Read Groups (required by GATK)");
                    gatkInputs[i] = bamFile.getProvisionedPath();
                    Job bamIndexJob = workflow.createBashJob("index_bams_" + i);
                    bamIndexJob.addFile(bamFile);
                    
                    /*bamIndexJob.setCommand(getWorkflowBaseDir() + "/bin/" + this.bundledJRE + "/bin/java"
                            + " -Xmx4000M"
                            + " -jar " + getWorkflowBaseDir() + "/bin/picard-tools-" + this.picardVersion + "/BuildBamIndex.jar"
                            + " INPUT=" + bamFile.getProvisionedPath());
                    */
                    bamIndexJob.setCommand(getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion
                                         + "/samtools index "
                                         + bamFile.getProvisionedPath());
                    bamIndexJob.setMaxMemory("4000");
                    fixerJobs.add(bamIndexJob);
                }

                Job jobGATK = workflow.createBashJob("call_snps_" + i);
                StringBuilder gatkCommand = new StringBuilder();

                gatkCommand.append(gatkJava)
                        .append(" -Xmx2g -Djava.io.tmpdir=").append(this.gatkDirs[i])
                        .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                        .append("-R ").append(this.genomeFile).append(" ")
                        .append("-T UnifiedGenotyper -I ")
                        .append(gatkInputs[i]).append(" ")
                        .append("-o ").append(this.vcfFiles[i]).append(" ")
                        .append("-filterRNC ")
                        .append("-stand_call_conf ").append(this.standCallConf).append(" ")
                        .append("-stand_emit_conf ").append(this.standEmitConf).append(" ")
                        .append("-dcov ").append(this.dcov).append(" ")
                        .append("-L ").append(this.checkedSNPs);
                newVcfs++; // Used to track the number of newly generated vcf files and
                // provisioning .vcf.gz and .vcf.tbi files

                jobGATK.setCommand(gatkCommand.toString());
                jobGATK.setMaxMemory(getProperty("gatk_memory"));

                if (!this.queue.isEmpty()) {
                    jobGATK.setQueue(this.queue);
                }

                for (Job fixerJob : fixerJobs) {
                    jobGATK.addParent(fixerJob);
                }
                gatkJobs.add(jobGATK);

                Job jobGATK2 = workflow.createBashJob("calculate_depth_" + i);
                StringBuilder gatk2Command = new StringBuilder();

                gatk2Command.append(gatkJava).append(" -Xmx3g -Djava.io.tmpdir=").append(this.gatkDirs[i])
                        .append(" -jar ").append(getWorkflowBaseDir()).append("/bin/GenomeAnalysisTK-").append(this.gatkVersion).append("/GenomeAnalysisTK.jar ")
                        .append("-R ").append(this.genomeFile).append(" ")
                        .append("-T DepthOfCoverage ")
                        .append("-I ").append(gatkInputs[i]).append(" ")
                        .append("-o ").append(this.tempDir).append(this.makeBasename(this.bamFiles[i])).append(" ")
                        .append("-filterRNC ")
                        .append("-L ").append(this.checkedSNPs);

                jobGATK2.setCommand(gatk2Command.toString());
                jobGATK2.setMaxMemory(getProperty("gatk_memory"));

                if (!this.queue.isEmpty()) {
                    jobGATK2.setQueue(this.queue);
                }

                for (Job fixerJob : fixerJobs) {
                    jobGATK2.addParent(fixerJob);
                }
                //Additional step to create depth of coverage data - for individual fingerprint image generation
                Job jobFin = workflow.createBashJob("make_fingerprint_file_" + i);
                StringBuilder finCommand = new StringBuilder();

                String basename = this.baseNames.get(i);
                finCommand.append(getWorkflowBaseDir()).append("/dependencies/create_fin.pl")
                        .append(" --refvcf=").append(this.checkedSNPs)
                        .append(" --genotype=").append(this.vcfFiles[i])
                        .append(" --coverage=").append(this.tempDir).append(basename).append(".sample_interval_summary")
                        .append(" --datadir=").append(this.tempDir)
                        .append(" --outdir=").append(this.dataDir).append(this.finDir)
                        .append(" --basename=").append(basename);
                // provision .fin file
                jobFin.addFile(this.createOutputFile(this.dataDir + this.finDir + basename + ".fin", "text/plain", this.manualOutput));

                jobFin.setCommand(finCommand.toString());
                jobFin.setMaxMemory("4000");
                jobFin.addParent(jobGATK2);

                if (!this.queue.isEmpty()) {
                    jobFin.setQueue(this.queue);
                }

                gatkJobs.add(jobFin);
            }

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

            for (String base : baseNames) {
                jobVcfPrep.addFile(this.createOutputFile(this.dataDir + base + SNPFILE_SUFFIX + ".gz.tbi", "application/tbi", this.manualOutput));
                jobVcfPrep.addFile(this.createOutputFile(this.dataDir + base + SNPFILE_SUFFIX + ".gz", "application/vcf-4-gzip", this.manualOutput));
            }

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
