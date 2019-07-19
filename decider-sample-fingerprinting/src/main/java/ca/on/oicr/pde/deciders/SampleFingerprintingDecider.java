/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 * Decider modified to work with Sample Fingerprinting workflow 2.0
 *
 * @author pruzanov@oicr.on.ca
 */
public class SampleFingerprintingDecider extends OicrDecider {

    private final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private String studyName = "";
    private String watchersList = "";
    private String customVcfList = "";
    private String queue = " ";

    //these params should come from settings xml file
    private String genomeFile;

    //Previous workflow runs and input files
    private String existingMatrix = "";
    private String templateTypeFilter = "";
    private String reseqTypeFilter = "";
    private boolean separate_platforms = true;
    private boolean allow_singletons   = false;
    private boolean mixedCoverageMode  = false;
    private Map<String, Map> reseqType;
    private String SNPConfigFile;
    private Map<String, BeSmall> fileSwaToSmall;

    //Static strings
    private final static String VCF_EXT = ".snps.raw.vcf.gz";
    private final static String TBI_EXT = ".snps.raw.vcf.gz.tbi";
    private final static String FIN_EXT = ".fin";
    private static final String VCF_GZIP_METATYPE = "application/vcf-4-gzip";
    private static final String TBI_METATYPE      = "application/tbi";
    private static final String FIN_METATYPE      = "text/plain";
    private static final String DEFAULT_SNPCONFIG_FILE = "/.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml";
    private static final String HOTSPOTS_DB_TOKEN = "file.hotspots_file";
    
    public SampleFingerprintingDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        defineArgument("ini-file", "Optional: the location of the INI file.", false);
        defineArgument("template-type", "Optional: name of the study that we need to analyze.", false);
        defineArgument("resequencing-type", "Optional: resequencing type for templates other than WG", false);
        defineArgument("existing-matrix", "Optional: existing matrix from previous workflow run(s)", false);
        defineArgument("queue", "Optional: Set the queue (for example, to production)", false);
        defineArgument("output-path", "Optional: the path where the files should be copied to "
                     + "after analysis. Corresponds to output-prefix in INI file. Default: ./", false);
        defineArgument("output-folder", "Optional: the name of the folder to put the output into relative to "
                     + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results", false);
        defineArgument("config-file", "Optional. Path to a config file in .xml format "
                     + "Default: /.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml", false);
        defineArgument("watchers-list", "Optional: Comma-separated list of emails for people interested in monitoring this workflow", false);
        defineArgument("allow-singletons", "Optional: A boolean flag that control inclusion singletons (donors with one bam file) in the final report", false);
        defineArgument("mixed-coverage", "Optional: A boolean flag that defines the algorithm of jaccard matrix calculation", false);
        defineArgument("manual-output", "Optional: Set the manual output. Default: false", false);
        defineArgument("separate-platforms", "Optional: Separate sequencing platforms, i.e. MiSeq and HiSeq. Default: true", false);
    }

    @Override
    public ReturnValue init() {

        Log.debug("INIT");
        String[] metaTypes = {TBI_METATYPE, VCF_GZIP_METATYPE, FIN_METATYPE};
        this.setMetaType(Arrays.asList(metaTypes));

        this.reseqType = new HashMap<String, Map>();
        // Group by template type if no other grouping selected
        if (!this.options.has("group-by")) {
            if (this.options.has("study-name")) {
                this.setGroupingStrategy(Header.STUDY_TITLE);
            }
            if (this.options.has("root-sample-name")) {
                // TODO: net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header needs to be updated to support ROOT_SAMPLE_SWA
                // uncomment when fixed... this.setGroupingStrategy(Header.ROOT_SAMPLE_SWA);
                // error out until the above is implemented
                throw new RuntimeException("ROOT_SAMPLE_SWA needs to be implemented in FindAllTheFiles.Header");
            }
            if (this.options.has("sequencer-run-name")) {
                this.setGroupingStrategy(Header.SEQUENCER_RUN_NAME);
            }
        } else {
            Log.warn("Passing group-by parameter overrides the defaults, I hope you know what you are doing");
        }

        if (this.options.has("mixed-coverage")) {
            this.mixedCoverageMode = Boolean.valueOf(options.valueOf("mixed-coverage").toString());
            Log.debug("Setting mixed coverage mode, default is false");
	}
        
        if (this.options.has("separate-platforms")) {
            String newSepValue = options.valueOf("separate-platforms").toString();
            if (newSepValue.equalsIgnoreCase("false")) {
               this.separate_platforms = false; 
            } else {
                Log.debug("Invalid setting for separate-platforms, using default [true]");
            }
	}
        
        if (this.options.has("allow-singletons")) {
            String newSepValue = options.valueOf("allow-singletons").toString();
            if (newSepValue.equalsIgnoreCase("true")) {
               this.allow_singletons = true; 
            } else {
                Log.debug("Invalid setting for allow-singletons, using default [false]");
            }
	}

        if (this.options.has("existing-matrix")) {
            this.existingMatrix = options.valueOf("existing-matrix").toString();
            if (null == this.existingMatrix || this.existingMatrix.isEmpty()) {
                this.existingMatrix = "";
            }
        }

        if (this.options.has("watchers-list")) {
            String commaSepWatchers = options.valueOf("watchers-list").toString();

            if (null != commaSepWatchers && !commaSepWatchers.isEmpty()) {
                Log.warn("We have " + commaSepWatchers + " for watchers");
                String[] watchers = commaSepWatchers.split(",");
                for (String email : watchers) {
                        this.watchersList += this.watchersList.isEmpty() ? email : "," + email;
                }
            }
            if (this.watchersList.isEmpty() || !this.watchersList.contains("@")) {
                this.watchersList = "";
            }
        }

        if (this.options.has("study-name")) {
            this.studyName = options.valueOf("study-name").toString();
            if (this.studyName.contains(" ")) {
                this.studyName = this.studyName.replaceAll(" ", "_");
            }
        } else {
            Log.warn("study-name parameter is not set, will try to determine it automatically");
        }

        if (this.options.has("config-file")) {
            this.SNPConfigFile = options.valueOf("config-file").toString();
            Log.warn("Got custom config file, will use settings from user's input");
        } else {
            this.SNPConfigFile = DEFAULT_SNPCONFIG_FILE;
            Log.warn("Using default config file " + DEFAULT_SNPCONFIG_FILE);
        }

        if (this.options.has("genome-file")) {
            this.genomeFile = options.valueOf("genome-file").toString();
            Log.warn("Got genome file, will use it instead of the dafault one (hg19)");
        }

        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
        } else {
            Log.debug("template-type parameter ids not set, will include all available types from the study");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        } else {
            Log.debug("queue parameter is not set, will not set this parameter");
        }
        
        // Setting resequencing type works like setting up a filter, the decider should fiure out reseq. type automatically
        if (this.options.has("resequencing-type")) {
            this.reseqTypeFilter = this.options.valueOf("resequencing-type").toString();
        } else {
            Log.debug("resequencing type parameter ids not set, will include all available types from the study");
        }

        // allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }
       
    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {

        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String targetTemplateType     = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // If nulls, set to NA
        if (null == targetResequencingType || targetResequencingType.isEmpty()) {
            targetResequencingType = "NA";
        }
        if (null == targetTemplateType || targetTemplateType.isEmpty()) {
            targetTemplateType = "NA";
        }
        
        // Get config if don't have it yet
        if (!this.reseqType.containsKey(targetTemplateType + targetResequencingType)) {
            boolean refsOK = this.configFromParsedXML(this.SNPConfigFile, targetTemplateType, targetResequencingType);
            if (!refsOK) {
                Log.error("References are not set for " + targetResequencingType + ", skipping");
                return false;
            }
        }
        
        //GP-470 Check hotspot file attribute - don't use the file if hotspots list is not matching the one we want
        String allowedHotspots = this.reseqType.get(targetTemplateType + targetResequencingType).get("file").toString();

        if (null != returnValue.getAttribute(HOTSPOTS_DB_TOKEN)) {
            if (!returnValue.getAttribute(HOTSPOTS_DB_TOKEN).equals(allowedHotspots)) {
                return false;
            }
        } else {
            // We can only assume that hotspots (if not annotated) will be from the
            // original set therefore we should fail all files without hotspots annotated
            // if current hotspots are not from the default file
            if (!this.SNPConfigFile.equals(DEFAULT_SNPCONFIG_FILE)) {
                return false;
            }
        }
                              
        // Check filters
        if (!this.reseqTypeFilter.isEmpty() && !this.reseqTypeFilter.equals(targetResequencingType)) {
            return false;
        }
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equals(targetTemplateType)) {
            return false;
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {

        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        //group files according to the designated header (e.g. sample SWID)
        List myTypes = this.getMetaType();

        //Make first pass, collect files in Set objects
        Set<String> tbiChecker = new HashSet<String>();
        Set<String> finChecker = new HashSet<String>();

        for (ReturnValue r : vals) {
            String fileMIME = r.getFiles().get(0).getMetaType();
            String filePath = r.getFiles().get(0).getFilePath();
            if (!myTypes.contains(fileMIME)) {
                continue;
            }
            if (fileMIME.equals(TBI_METATYPE)) {
                tbiChecker.add(this.makeBasePath(filePath, TBI_EXT));
            }
            if (fileMIME.equals(FIN_METATYPE)) {
                finChecker.add(this.makeBasePath(filePath, FIN_EXT));
            }
        }

        for (ReturnValue r : vals) {
            String currentRV = r.getAttributes().get(groupBy);
            String template_type = r.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            if (null == this.studyName || this.studyName.isEmpty()) {
                FileAttributes fa = new FileAttributes(r, r.getFiles().get(0));
                String t = fa.getDonor();
                this.studyName = t.substring(0, t.indexOf("_"));
                Log.debug("Extracted Study Name " + this.studyName);
            }
            if (null == currentRV || null == template_type || (!this.templateTypeFilter.isEmpty() && !template_type.equals(this.templateTypeFilter))) {
                continue;
            }

            //Make sure that BeSmall gets only vcf files
            if (!r.getFiles().get(0).getMetaType().equals(VCF_GZIP_METATYPE)) {
                continue;
            }

            String filePath = this.makeBasePath(r.getFiles().get(0).getFilePath(), VCF_EXT);
            //Make sure we have proper sets of data files
            if (!tbiChecker.contains(filePath) || !finChecker.contains(filePath)) {
                continue;
            }

            BeSmall currentSmall = new BeSmall(r);
            fileSwaToSmall.put(r.getAttribute(Header.FILE_SWA.getTitle()), currentSmall);

            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getPath());
                iusDeetsToRV.put(fileDeets, r);
            } 
            //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);

                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\t instead of file "
                            + "\n\t" + oldSmall.getDate());
                    iusDeetsToRV.put(fileDeets, r);
                } else {
                    Log.debug("Disregarding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\tas older than duplicate sequencer run/lane/barcode in favour of "
                            + "\n\t" + oldSmall.getDate());
                    Log.debug(currentDate + " is before " + oldDate);
                }
            }

        }

        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (in the case of this workflow, by template type)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).groupByAttribute;
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }
    
    @Override
    public ReturnValue customizeRun(WorkflowRun run) {
        Log.debug("INI FILE:" + run.getIniFile().toString());
        String checkedSNPs = "";
        String checkPoints = "";
        List<String> vcfList = new ArrayList<String>();
        
        for (FileAttributes atts : run.getFiles()) {
            if (checkedSNPs.isEmpty()) {
                String fileKey = fileSwaToSmall.get(atts.getOtherAttribute(Header.FILE_SWA.getTitle())).getTemplateType() +
                                 fileSwaToSmall.get(atts.getOtherAttribute(Header.FILE_SWA.getTitle())).getReseqType();
                if (this.reseqType.containsKey(fileKey)) {
                    checkedSNPs = this.reseqType.get(fileKey).get("file").toString();
                    checkPoints = this.reseqType.get(fileKey).get("points").toString();
                }
            }
         if (atts.getMetatype().equalsIgnoreCase(VCF_GZIP_METATYPE)) {
             vcfList.add(atts.getPath());
         }
        }

        StringBuilder inputString = new StringBuilder();
        // Check for name unique-ness, skip duplicates
        Set<String> uniqChecker = new HashSet<String>();
        String[] inputs = vcfList.toArray(new String[0]);
        for (String inputVcf : inputs) {
            if (uniqChecker.contains(inputVcf)) {
                Log.stderr("Found duplicate file in the list of extra inputs, will skip [" + inputVcf + "]");
                continue;
            }

            uniqChecker.add(inputVcf);
            if (!inputString.toString().isEmpty()) {
                inputString.append(",");
            }
            inputString.append(inputVcf);

        }

        if (!this.customVcfList.isEmpty()) {
            String[] inputsExtra = this.customVcfList.split(",");
            for (String extraInput : inputsExtra) {
                if (uniqChecker.contains(extraInput)) {
                    Log.stderr("Found duplicate file in the list of extra inputs, will skip [" + extraInput + "]");
                    continue;
                }
                inputString.append(",").append(extraInput);
            }
        }

        run.addProperty("input_files", inputString.toString());
        run.addProperty("mixed_coverage", Boolean.toString(this.mixedCoverageMode));
        run.addProperty("allow_singletons", Boolean.toString(this.allow_singletons));
        
        if (!this.queue.isEmpty()) {
          run.addProperty("queue", this.queue);
        } else {
          run.addProperty("queue", " ");
        }
        
        if (null != checkedSNPs && !checkedSNPs.isEmpty()) {
            run.addProperty("checked_snps", checkedSNPs);
            run.addProperty("check_points", checkPoints);
        }

        if (null != this.genomeFile) {
            run.addProperty("genome_file", this.genomeFile);
        }

        run.addProperty("study_name", this.studyName);

        if (!this.existingMatrix.isEmpty()) {
            run.addProperty("existing_matrix", this.existingMatrix);
        } else {
            run.addProperty("existing_matrix", " ");
        }

        if (!this.watchersList.isEmpty()) {
            run.addProperty("watchers_list", this.watchersList);
        }

        return new ReturnValue();
    }

    private boolean configFromParsedXML(String fileName, String templateType, String resequencingType) {

        try {
            File fXmlFile = new File(fileName);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);

	    //optional, but recommended
            //read this - http://stackoverflow.com/questions/13786607/normalization-in-dom-parsing-with-java-how-does-it-work
            doc.getDocumentElement().normalize();
            NodeList nList = doc.getElementsByTagName("template_type");
            for (int temp = 0; temp < nList.getLength(); temp++) {

                Node nNode = nList.item(temp);
                if (nNode.getNodeType() == Node.ELEMENT_NODE) {
                    Element nElement = (Element) nNode;
                    if (!templateType.equals(nElement.getAttribute("id"))) {
                        continue;
                    }
                }

                if (nNode.hasChildNodes()) {
                    NodeList children = nNode.getChildNodes();
                    for (int cld = 0; cld < children.getLength(); cld++) {
                        Node cNode = children.item(cld);
                        if (cNode.getNodeType() == Node.ELEMENT_NODE) {
                            Element cElement = (Element) cNode;
                            if (!resequencingType.equals(cElement.getAttribute("id"))) {
                                continue;
                            }
                            Map<String, String> rtypeData = new HashMap<String, String>();
                            rtypeData.put("file", cElement.getElementsByTagName("checked_snps").item(0).getTextContent());
                            rtypeData.put("points", cElement.getElementsByTagName("check_points").item(0).getTextContent());
                            this.reseqType.put(templateType.concat(resequencingType), rtypeData);
                            return true;
                        }
                    }
                }
            }
        } catch (FileNotFoundException fnf) {
            Log.error("File is not found");
        } catch (NullPointerException np) {
            Log.error("A value in config file for " + resequencingType + " is not initialized");
        } catch (NumberFormatException nf) {
            Log.error("A number of checkponts in the config file is not set properly, should be an integer");
        } catch (Exception e) {
            e.printStackTrace();
        }
        return false;
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(SampleFingerprintingDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String path = null;
        private String templateType;
        private String reseqType;
        
        public BeSmall(ReturnValue rv) {
            try {
                date = DATE_FORMAT.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            
            // Having metatype as part of details is needed since we deal with multiple mime types her
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + fa.getMetatype();
            // We are going to group by template type only (if we did not receive template type as a parameter)
            this.templateType = fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
            this.reseqType    = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
            
            StringBuilder groupBy = new StringBuilder(this.templateType);
             String platformID = rv.getAttribute("Sequencer Run Platform ID");
            if (null != this.reseqType) {
                groupBy.append(":").append(this.reseqType);
            } else {
                this.reseqType = "NA";
                groupBy.append(":").append(this.reseqType);
            }
            //Depending on user's options, platformID may be used for grouping
            if (separate_platforms)
                groupBy.append(":").append(platformID);
            this.setGroupByAttribute(groupBy.toString());
            
            path = rv.getFiles().get(0).getFilePath();
        }

        public Date getDate() {
            return date;
        }

        public final void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public final void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public final void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public final void setPath(String path) {
            this.path = path;
        }

        /**
         * @return the templateType
         */
        public String getTemplateType() {
            return templateType;
        }

        /**
         * @return the reseqType
         */
        public String getReseqType() {
            return reseqType;
        }
    }

    // Service Functions
    private String makeBasePath(String name, String ext) {
        return name.contains(ext) ? name.substring(0, name.lastIndexOf(ext)) : name;
    }

}
