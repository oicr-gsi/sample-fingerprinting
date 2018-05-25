/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;

/**
 * SampleFingerprinting v 2.0 - split from original (1.1) version and working
 * with Fingerprint Collector workflow
 *
 * @author pruzanov@oicr.on.ca
 */
public class SampleFingerprintingWorkflow extends OicrWorkflow {

    private String[] vcfFiles;
    private String existingMatrix = "";
    private String studyName = "";
    private String watchersList = "";
    private String dataDir;
    private String tabixVersion;
    private String vcftoolsVersion;
    private String checkPoints;
    private String checkedSnps;
    private String queue;
    private boolean manualOutput;
    private boolean allowSingletons;
    private boolean mixedCoverageMode = false;
    private int jChunkSize; // Optimal (for speed) allowed number of vcf files when jaccard_indexing step doesn't fork into multiple sub-jobs
    //Static integers
    private final int MIN_CHUNK_SIZE = 50;
    private final int MAX_CHUNKS     = 40;   // n in (n^2 + n)/2 which defines N of chunks (it should not be more than 900
    private final int MAX_FILES      = 500;  // Maximum number of files to pass to the script for sym.linking
    private final int MAX_INPUTS     = 3500; // Recommended max size of the data set
    //Static strings
    private final String finDir           = "finfiles/";
    private final String reportName       = "sample_fingerprint";
    private final static String VCF_EXT   = ".snps.raw.vcf.gz";
    private final static String TABIX_EXT = ".snps.raw.vcf.gz.tbi";
    private final static String FIN_EXT   = ".fin";
    private static final String HOTSPOTS_TOKEN = "hotspots_file";

    @Override
    public Map<String, SqwFile> setupFiles() {
        // Set up reference, bam and vcf files here     
        try {
            if (getProperty("input_files") == null) {
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.SEVERE, "input_files is not set, we need at least one bam file");
                return (null);
            } else {
                this.vcfFiles = getProperty("input_files").split(",");
                if (this.vcfFiles.length > MAX_INPUTS)
                    Log.stderr("Number of inputs files exceeds recommended threshold (" + MAX_INPUTS + ") which may be impractical for this workflow");
                this.jChunkSize = MIN_CHUNK_SIZE;
                int chunks = this.vcfFiles.length / this.jChunkSize;
                if (chunks > MAX_CHUNKS) {
                    this.jChunkSize = this.vcfFiles.length/MAX_CHUNKS;
                }
            }

            this.checkPoints     = getProperty("check_points");           
            this.existingMatrix  = getOptionalProperty("existing_matrix","");
            this.watchersList    = getOptionalProperty("watchers_list", "");
            this.tabixVersion    = getProperty("tabix_version");
            this.vcftoolsVersion = getProperty("vcftools_version");
            this.studyName       = getOptionalProperty("study_name", "");
            this.queue           = getOptionalProperty("queue", "");
            this.checkedSnps     = getProperty("checked_snps");

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
            
            if (getProperty("mixed_coverage") != null ) {
               try {
                   this.mixedCoverageMode = Boolean.valueOf(getProperty("mixed_coverage"));
               } catch (NumberFormatException ne) {
                   this.mixedCoverageMode = false;
               }
            }

            if (getProperty("allow_singletons") == null) {
                this.allowSingletons = false;
                Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, "Setting allow singletons to default (false)");
            } else {
                String manualCheck = getOptionalProperty("allow_singletons", "false");
                try {
                    this.allowSingletons = Boolean.valueOf(manualCheck);
                } catch (NumberFormatException e) {
                    this.allowSingletons = false;
                }
            }

            // Set up all inputs, assume that all vcf.gz, vcf.gz.tbi and .fin files sit in triplets in the same directories
            for (int i = 0; i < this.vcfFiles.length; i++) {
                //Using first file, try to guess the study name if it was not provided as an argument in .ini file
                if (i == 0 && this.studyName.isEmpty()) {
                    if (vcfFiles[i].matches("SWID_\\d+_\\D+_\\d+")) {
                        String tempName = vcfFiles[i].substring(vcfFiles[i].lastIndexOf("SWID_"));
                        this.studyName = tempName.substring(0, tempName.indexOf("_") - 1);
                        break;
                    } else {
                        this.studyName = "UNDEF";
                    }
                }

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
            //Setup data and finfiles dirs
            this.dataDir = getOptionalProperty("data_dir", "data/");
            if (!this.dataDir.endsWith("/")) {
                this.dataDir += "/";
            }
            this.addDirectory(this.dataDir);
            this.addDirectory(this.dataDir + this.finDir);
        } catch (Exception e) {
            Logger.getLogger(SampleFingerprintingWorkflow.class.getName()).log(Level.WARNING, null, e);
        }

    }

    @Override
    public void buildWorkflow() {
        try {
            //Need the decider to check the headers of the .bam files for the presence of RG and abscence of empty Fields (Field: )
            Workflow workflow = this.getWorkflow();
            List<Job> upstreamJobs = new ArrayList<Job>();
            List<Job> linkerJobs   = new ArrayList<Job>();
            // A small copy job
            Job job_copy = workflow.createBashJob("copy_resources");
            job_copy.setCommand("cp -R " + getWorkflowBaseDir()
                    + "/data/* "
                    + this.dataDir);

            job_copy.setMaxMemory(getProperty("memory_job_copy"));
            if (!this.queue.isEmpty()) {
                job_copy.setQueue(this.queue);
            }
            upstreamJobs.add(job_copy);

            final String extList = VCF_EXT + "," + FIN_EXT + "," + TABIX_EXT;
            for (int l = 0; l < this.vcfFiles.length; l += MAX_FILES) {
                int endOfLoop = l + MAX_FILES > this.vcfFiles.length ? this.vcfFiles.length : l + MAX_FILES - 1;
                StringBuilder basePathList = new StringBuilder(this.makeBasepath(this.vcfFiles[l], VCF_EXT));
                for (int v = l + 1; v < endOfLoop; v++) {
                    basePathList.append(",").append(this.makeBasepath(this.vcfFiles[v], VCF_EXT));

                }
                // Prepare the inputs (make symlinks in dataDir)
                Job job_link_maker = workflow.createBashJob("link_inputs");
                job_link_maker.setCommand(getWorkflowBaseDir() + "/dependencies/inputs_linker.pl"
                        + " --list "    + basePathList.toString()
                        + " --datadir " + this.dataDir
                        + " --findir "  + this.dataDir + this.finDir
                        + " --extensions " + extList);
                job_link_maker.setMaxMemory(getProperty("memory_job_link_maker"));
                if (!this.queue.isEmpty()) {
                    job_link_maker.setQueue(this.queue);
                }
                linkerJobs.add(job_link_maker);
            }

            String chunkList = null;
            
            if (!this.mixedCoverageMode) {
            if (this.existingMatrix.isEmpty() && this.vcfFiles.length > this.jChunkSize) {

                StringBuilder chunkedResults = new StringBuilder();
                int vcf_chunks = (int) Math.ceil((double) this.vcfFiles.length / (double) this.jChunkSize);
                int resultID = 1;

                // Combine chunks and run the jobs, registering results for later use
                for (int c = 0; c < vcf_chunks; c++) {
                    for (int cc = c; cc < vcf_chunks; cc++) {

                        Job job_list_writer = workflow.createBashJob("make_list" + resultID);
                        chunkList = "jaccard.chunk." + resultID + ".list";
                        job_list_writer.setCommand(getWorkflowBaseDir() + "/dependencies/write_list.pl"
                                + " --datadir "    + this.dataDir
                                + " --segments \"" + c * this.jChunkSize + ":" + ((c + 1) * this.jChunkSize - 1) + ","
                                + cc * this.jChunkSize + ":" + ((cc + 1) * this.jChunkSize - 1) + "\""
                                + " > " + this.dataDir + chunkList);

                        job_list_writer.setMaxMemory(getProperty("memory_job_list_writer"));
                        if (!this.queue.isEmpty()) {
                            job_list_writer.setQueue(this.queue);
                        }

                        for (Job j : linkerJobs) {
                            job_list_writer.addParent(j);
                        }

                        Job job_jchunk = workflow.createBashJob("make_matrix_" + resultID);
                        String chunkName = "jaccard.chunk." + resultID + ".csv";
                        String sep = chunkedResults.toString().isEmpty() ? "" : ",";
                        chunkedResults.append(sep).append(this.dataDir).append(chunkName);

                        resultID++;

                        job_jchunk.setCommand(getWorkflowBaseDir() + "/dependencies/jaccard_coeff.matrix.pl"
                                + " --list "        + this.dataDir + chunkList
                                + " --vcf-compare " + getWorkflowBaseDir() + "/bin/vcftools_" + this.vcftoolsVersion + "/bin/vcf-compare"
                                + " --datadir "     + this.dataDir
                                + " --tabix "       + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion
                                + " --studyname "   + this.studyName
                                + " > " + this.dataDir + chunkName);

                        job_jchunk.setMaxMemory(getProperty("memory_job_jchunk"));
                        if (!this.queue.isEmpty()) {
                            job_jchunk.setQueue(this.queue);
                        }

                        job_jchunk.addParent(job_list_writer);
                        upstreamJobs.add(job_jchunk);
                    }
                }

                // The result of the last job will be existingMatrix
                this.existingMatrix = chunkedResults.toString();
            }

            // Next job is jaccard_coeff.matrix script, need to create|update giant matrix of jaccard indexes
            // using vcftools, colors should not be assigned at his step. Will be final jaccard index job          
            Job job_list_writer2 = workflow.createBashJob("make_Final_list");
            chunkList = "jaccard.chunks.list";
            job_list_writer2.setCommand(getWorkflowBaseDir() + "/dependencies/write_list.pl"
                    + " --datadir " + this.dataDir
                    + " > "         + this.dataDir + chunkList);

            job_list_writer2.setMaxMemory(getProperty("memory_job_list_writer2"));
            if (!this.queue.isEmpty()) {
                job_list_writer2.setQueue(this.queue);
            }

            for (Job upstreamJob : upstreamJobs) {
                job_list_writer2.addParent(upstreamJob);
            }
            upstreamJobs.clear();
            upstreamJobs.add(job_list_writer2);
        }

            //Similarity-calculating job, operates on chunks or list of finfiles
            SqwFile matrix = this.createOutputFile(this.dataDir + this.studyName + "_jaccard.matrix.csv", "text/plain", this.manualOutput);
            Job job_jaccard = workflow.createBashJob("make_matrix");           

            if (this.mixedCoverageMode) {
                job_jaccard.setCommand(getWorkflowBaseDir() + "/dependencies/jaccard_coeff.matrix.mc.pl"
                        + " --findir "  + this.dataDir + this.finDir 
                        + " > " + matrix.getSourcePath());
            } else {
                job_jaccard.setCommand(getWorkflowBaseDir() + "/dependencies/jaccard_coeff.matrix.pl"
                        + " --list "        + this.dataDir + chunkList
                        + " --vcf-compare " + getWorkflowBaseDir() + "/bin/vcftools_" + this.vcftoolsVersion + "/bin/vcf-compare"
                        + " --datadir "     + this.dataDir
                        + " --tabix "       + getWorkflowBaseDir() + "/bin/tabix-" + this.tabixVersion
                        + " --studyname "   + this.studyName
                        + " > " + matrix.getSourcePath());
                if (!this.existingMatrix.isEmpty()) {
                    job_jaccard.getCommand().addArgument("--existing_matrix " + this.existingMatrix);
                }
            }

            job_jaccard.setMaxMemory(getProperty("memory_job_jaccard"));
            if (!this.queue.isEmpty()) {
                job_jaccard.setQueue(this.queue);
            }
            matrix.getAnnotations().put(HOTSPOTS_TOKEN, this.checkedSnps);
            job_jaccard.addFile(matrix);
            for (Job upstreamJob : upstreamJobs) {
                job_jaccard.addParent(upstreamJob);
            }

            // Images generated here: matrix slicing/color assignment, flagging suspicious files, wrapper for R script
            Job make_pics = workflow.createBashJob("make_report");
            StringBuilder reportCommand = new StringBuilder();
            reportCommand.append("perl ").append(getWorkflowBaseDir()).append("/dependencies/make_report.pl")
                    .append(" --matrix ").append(matrix.getSourcePath())
                    .append(" --refsnps ").append(this.checkPoints)
                    .append(" --tempdir ").append(this.dataDir).append(this.finDir)
                    .append(" --datadir ").append(this.dataDir)
                    .append(" --studyname ").append(this.studyName);
            
            if (this.allowSingletons) {reportCommand.append(" --allow-singletons "); }               
            reportCommand.append(" > ").append(this.dataDir).append("index.html");

            make_pics.setCommand(reportCommand.toString());
            make_pics.addParent(job_jaccard);
            make_pics.setMaxMemory(getProperty("memory_make_pics"));
            if (!this.queue.isEmpty()) {
                make_pics.setQueue(this.queue);
            }

            // Make report table here using index.html and jaccard_coeff.matrix.csv
            Job prox_table = workflow.createBashJob("proximity_table");
            prox_table.setCommand("perl " + getWorkflowBaseDir() + "/dependencies/make_table.pl"
                                + " --index "  + this.dataDir + "index.html"
                                + " --matrix " + matrix.getSourcePath()
                                + " --out " + this.dataDir + this.studyName + ".proximity_table.csv");
            if (this.allowSingletons) {prox_table.getCommand().addArgument(" --singletons "); }

            prox_table.addParent(make_pics);
            prox_table.setMaxMemory(getProperty("memory_prox_table"));
            if (!this.queue.isEmpty()) {
                prox_table.setQueue(this.queue);
            }

            // Zip finfiles and similarity matrix for customization in the webtool
            Job zip_fins = workflow.createBashJob("zip_finfiles");
            zip_fins.setCommand("zip -r " + this.dataDir + "customize.me.zip "
                    + this.dataDir + this.finDir + " "
                    + matrix.getSourcePath());
            zip_fins.addParent(prox_table);
            zip_fins.setMaxMemory(getProperty("memory_zip_fins"));
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
                    + this.dataDir + this.studyName + ".proximity_table.csv "
                    + this.dataDir + "customize.me.zip "
                    + this.dataDir + "*.html");
            zip_report.addParent(zip_fins);
            zip_report.setMaxMemory(getProperty("memory_zip_report"));

            report_file.getAnnotations().put(HOTSPOTS_TOKEN, this.checkedSnps);
            zip_report.addFile(report_file);

            if (!this.queue.isEmpty()) {
                zip_report.setQueue(this.queue);
            }

            // Next job is report-emailing script (check html for flagged samples, send email(s) to interested people)
            if (!this.watchersList.isEmpty()) {
                Job job_alert = workflow.createBashJob("send_alert");
                job_alert.setCommand("perl " + getWorkflowBaseDir() + "/dependencies/plotReporter.pl"
                        + " --watchers "  + this.watchersList
                        + " --report "    + this.dataDir + "index.html"
                        + " --studyname " + this.studyName
                        + " --bundle "    + report_file.getOutputPath());
                job_alert.setMaxMemory(getProperty("memory_job_alert"));
                job_alert.addParent(zip_report);
                if (!this.queue.isEmpty()) {
                    job_alert.setQueue(this.queue);
                }
            }

        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }

    }

    private String makeBasepath(String path, String extension) {
        return path.substring(0, path.lastIndexOf(extension));
    }
}
