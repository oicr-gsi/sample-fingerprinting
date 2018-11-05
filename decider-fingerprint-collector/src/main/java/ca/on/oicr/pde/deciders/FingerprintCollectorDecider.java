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
import java.util.List;
import java.util.Map;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.common.util.maptools.MapTools;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

public class FingerprintCollectorDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;
    private String gatkMemory;

    private String groupBy = "";
    private String customBamList = "";
    //GATK:
    private String gatkPrefix = "./";
    private double standCallConf;
    private double standEmitConf;
    private int dcov;

    private String queue = " ";
    private final static String BAM_METATYPE = "application/bam";
    private final static int HUMAN_ORG_ID = 31;
   
    private String templateTypeFilter   = "";
    private String resequenceTypeFilter = "";
    private String studyName;
    private boolean preprocessBam = false;
    private Map<String, Map> reseqType;
    private String SNPConfigFile;
    private static final String DEFAULT_SNPCONFIG_FILE = "/.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml";

    public FingerprintCollectorDecider () {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        defineArgument("preprocess-bam", "Optional. Set the flag that tells the workflow to run bam re-ordering "
                + "and adding RG (read groups) which may be absent in some cases", false);
        defineArgument("template-type", "Optional. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type", false);
        defineArgument("resequencing-type", "Optional. Set the resequencing type filter to limit the workflow run "
                + "so that it runs on data only of this resequencing type", false);
        defineArgument("config-file", "Optional. Path to a config file in .xml format "
                     + "Default: /.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml", false);
        defineArgument("queue", "Optional: Set the queue (for example, to production)", false);
        defineArgument("gatk-memory", "Optional. Set the memory allocated to GATK jobs "
                + "when running the workflow, the default is 5000", false);
        defineArgument("stand-call-conf", "Optional. Set GATK parameter stand_call_conf "
                + "when running the workflow, the default is 50.0", false);
        defineArgument("stand-emit-conf", "Optional. Set GATK parameter stand_emit_conf "
                + "when running the workflow, the default is 10.0", false);
        defineArgument("dcov", "Optional. Set GATK parameter dcov when running the workflow, the default is 200", false);
        defineArgument("gatk-prefix", "Optional: the path to a dir on a low-latency filesystem for writing " 
                + "GATK temporary data. May prevent possible failures of a workflow run. Default: ./", false);

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(BAM_METATYPE));
        this.reseqType = new HashMap<String, Map>();
        //Handle .ini file - we accept only memory size allocated to different steps
        if (options.has("ini-file")) {
            File file = new File(options.valueOf("ini-file").toString());
            if (file.exists()) {
                String iniFile = file.getAbsolutePath();
                Map<String, String> iniFileMap = new HashMap<String, String>();
                MapTools.ini2Map(iniFile, iniFileMap);

                this.standCallConf = Double.valueOf(iniFileMap.get("stand_call_conf"));
                this.standEmitConf = Double.valueOf(iniFileMap.get("stand_emit_conf"));
                this.dcov          = Integer.valueOf(iniFileMap.get("dcov"));
                this.customBamList = iniFileMap.get("input_files");
                if (!this.customBamList.isEmpty())
                    Log.warn("Note that list of .bam files from this .ini will be APPENDED to the listt retrieved from metaDB");
            } else {
                Log.error("The given INI file does not exist: " + file.getAbsolutePath());
                System.exit(1);
            }

        }

        //Group by
        if (this.options.has("group-by")) {
            this.groupBy = options.valueOf("group-by").toString();
            if (!this.groupBy.equals(Header.FILE_SWA.toString()) &&
                !this.groupBy.equals(Header.SAMPLE_NAME.toString()) &&
                !this.groupBy.equals(Header.SEQUENCER_RUN_NAME.toString())) {
                Log.debug("Unsupported group-by parameter passed, will revert to default group-by FILE_SWA");
                this.groupBy = Header.FILE_SWA.toString();
            } else {
              if (!this.groupBy.equals(Header.FILE_SWA.toString())) {
                  Log.debug("group-by [" + this.groupBy + "] passed, it will override the default (FILE_SWA)");
              }
            }
        }
        
        //Very important - picking config file for hotspot locations
        if (this.options.has("config-file")) {
            this.SNPConfigFile = options.valueOf("config-file").toString();
        } else {
            this.SNPConfigFile = DEFAULT_SNPCONFIG_FILE;
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
            Log.debug("Setting template type is not necessary, however if set the decider will run the workflow only on this type of data");
        }

        if (this.options.has("resequencing-type")) {
            this.resequenceTypeFilter = options.valueOf("resequencing-type").toString();
            Log.debug("Setting resequencing type is not necessary, however if set the decider will run the workflow only on this type of data");
        }
        
        if (this.options.has("gatk-memory")) {

            this.gatkMemory = options.valueOf("gatk-memory").toString();
            Log.debug("Setting memory allocated for GATK step to " + this.gatkMemory + " Megabytes as requested");
        }

        if (this.options.has("gatk-prefix")) {
          this.gatkPrefix = options.valueOf("gatk-prefix").toString();
          if (this.gatkPrefix.isEmpty())
	      this.gatkPrefix = "./";
          else if (!this.gatkPrefix.endsWith("/"))
              this.gatkPrefix += "/";
        }
        
        if (this.options.has("stand-call-conf")) {
            try {
                this.standCallConf = Double.valueOf(options.valueOf("stand-conf-call").toString());
            } catch (NumberFormatException nf) {
                Log.error("stand-conf-call failed to pass as double, make sure you supply a valid number");
                System.exit(1);
            }
        }

        if (this.options.has("stand-emit-conf")) {
            try {
                this.standEmitConf = Double.valueOf(options.valueOf("stand-emit-call").toString());
            } catch (NumberFormatException nf) {
                Log.error("stand-emit-call failed to pass as double, make sure you supply a valid number");
                System.exit(1);
            }
        }

        if (this.options.has("dcov")) {
            try {
                this.dcov = Integer.valueOf(options.valueOf("dcov").toString());
            } catch (NumberFormatException nf) {
                Log.error("dcov failed to pass as int, make sure you supply a valid number");
                System.exit(1);
            }
        }
        
        if (this.options.has("study-name")) {
            this.studyName = options.valueOf("study-name").toString();
	} else {
            Log.warn("study-name parameter is not set, will try to determine it automatically");
        }
        
        if (this.options.has("preprocess-bam")) {
            this.preprocessBam = Boolean.valueOf(options.valueOf("preprocess-bam").toString());
        }

        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }
    
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveRStype = false;

        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                String fileKey = bs.getReseqTemplateID();
                if (this.reseqType.containsKey(fileKey)) {
                    haveRStype = true;
                } else {
                    haveRStype = false;
                    break; // even if one of the files does not have reseq type supported, break
                }
            }
        }

        if (haveRStype) {
            return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        }

        Log.error("Resequencing Data for at least some of the inputs are not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String targetTemplateType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // If nulls, set to NA
        if (null == targetTemplateType || targetTemplateType.isEmpty()) {
            targetTemplateType = "NA";
        }
        if (null == targetResequencingType || targetResequencingType.isEmpty()) {
            targetResequencingType = "NA";
        }

        if (!fm.getMetaType().equals(BAM_METATYPE)) {
                return false;
        }
        // Check filters
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equals(targetTemplateType)) {
            return false;
        }
        if (!this.resequenceTypeFilter.isEmpty() && !this.resequenceTypeFilter.equals(targetResequencingType)) {
            return false;
        }
        //Check organism GP-473 check if it is indeed an integer
        String organismId  = returnValue.getAttribute("Sample Organism ID");
        if (null != organismId && !organismId.isEmpty()) {
            try {
                int passedId = Integer.valueOf(organismId);
                if (passedId != HUMAN_ORG_ID) {
                    Log.error("Organism other than H.sapience is not supported");
                    return false;
                }
            } catch (NumberFormatException nfe) {
                Log.error("Invalid value passed as Organism ID, assuming H.sapience");
            }   
        }

        // Get config if don't have it yet
        if (!this.reseqType.containsKey(targetTemplateType + targetResequencingType)) {
            boolean refsOK = this.configFromParsedXML(this.SNPConfigFile, targetTemplateType, targetResequencingType);
            if (!refsOK) {
                Log.error("References are not set for " + targetResequencingType + ", skipping");
                return false;
            }
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        //group files according to the designated header (e.g. sample SWID)
        List myTypes = this.getMetaType();
        for (ReturnValue r : vals) {
            if (!myTypes.contains(r.getFiles().get(0).getMetaType()))
                continue;

            String currentRV  = r.getAttributes().get(groupBy);
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
                

           // String groupByThis = this.groupBy.isEmpty() ? "" : this.handleGroupByAttribute(this.groupBy);
            BeSmall currentSmall = new BeSmall(r, this.groupBy);
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
            
            String currVal;
            if (groupBy.equals(Header.FILE_SWA.getTitle())) {
              currVal = r.getAttribute(Header.FILE_SWA.getTitle());
            } else {
              currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).groupByAttribute;
            }
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }
        
     return map;
    }

    @Override //Not used here
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        return a;
    }

    @Override
    public ReturnValue customizeRun(WorkflowRun run) {
       
        StringBuilder inputFiles  = new StringBuilder();
        StringBuilder rgDetails  = new StringBuilder();
        String checkedSNPs = "";
        String checkPoints = "";
        
        for (FileAttributes atts : run.getFiles()) {
          if (checkedSNPs.isEmpty()) {
           String fileKey = fileSwaToSmall.get(atts.getOtherAttribute(Header.FILE_SWA.getTitle())).getReseqTemplateID();
           if (this.reseqType.containsKey(fileKey)) {
             checkedSNPs = this.reseqType.get(fileKey).get("file").toString();
             checkPoints = this.reseqType.get(fileKey).get("points").toString();
           } else {
                
           }
          }
          
          if (atts.getMetatype().equals(BAM_METATYPE)) {
            if (inputFiles.length() != 0) {
                inputFiles.append(",");
                rgDetails.append(",");
            }
              inputFiles.append(atts.getPath());
              rgDetails.append(fileSwaToSmall.get(atts.getOtherAttribute(Header.FILE_SWA.getTitle())).getRGdata());
          }
        }       
        
        if (!this.customBamList.isEmpty()) {   
            inputFiles.append(",").append(this.customBamList);
        }
        run.addProperty("input_files", inputFiles.toString());
        run.addProperty("rg_details", rgDetails.toString());
        run.addProperty("gatk_prefix",this.gatkPrefix);
        run.addProperty("preprocess_bam", Boolean.toString(this.preprocessBam));
        run.addProperty("checked_snps", checkedSNPs);
        run.addProperty("check_points", checkPoints);

               
        if (this.standCallConf > 0) {
          run.addProperty("stand_call_conf", "" + this.standCallConf);
        }
        if (this.standEmitConf > 0) {
          run.addProperty("stand_emit_conf", "" + this.standEmitConf);
        }
        if (this.dcov > 0) {
          run.addProperty("dcov", "" + this.dcov);
        }
        if (!this.studyName.isEmpty()) {
          run.addProperty("study_name", this.studyName);
        }
        if (!this.queue.isEmpty()) {
          run.addProperty("queue", this.queue);
        } else {
          run.addProperty("queue", " ");
        }
        
        return new ReturnValue();
    }

    private boolean configFromParsedXML(String fileName, String templateType, String resequencingType) {
        // SHOULD BE OK, though fields that need to be defined also need to be customized
        try {
            File fXmlFile = new File(fileName);
            DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
            DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
            Document doc = dBuilder.parse(fXmlFile);

	// optional, but recommended
        // read this - http://stackoverflow.com/questions/13786607/normalization-in-dom-parsing-with-java-how-does-it-work
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
        params.add(FingerprintCollectorDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

        private Date   date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String path = null;
        private String reseqTemplateID = null;
        
        private String rgId = null;
        private String rgPl = null;
        private String rgPu = null;
        private String rgLb = null;
        private String rgSm = null;

        public BeSmall(ReturnValue rv, String groupBy) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            this.rgId = rv.getAttribute(Header.FILE_SWA.getTitle());
            this.rgPl = rv.getAttribute("Sequencer Run Platform Name");
            String rsType = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
            if (rsType == null || rsType.isEmpty() || rsType.equals(" ")) {
                rsType = "NA";
            }
            
            if (null != this.rgPl && !this.rgPl.isEmpty() && this.rgPl.contains(" ")) {
                this.rgPl = this.rgPl.substring(0, this.rgPl.indexOf(" "));
            }

            this.rgPu = fa.getSequencerRun() + "_" + fa.getLane() + "_" + fa.getBarcode();
            this.rgLb = fa.getLibrarySample();
            this.rgSm = fa.getDonor();
            
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            
            reseqTemplateID = fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE) + rsType;
            groupByAttribute = iusDetails;
            // If groupBy
            if (!groupBy.isEmpty()) {
                if (groupBy.equals("SAMPLE_NAME")) {
                    groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
                } else if (groupBy.equals("SAMPLE_SWA")) {
                    groupByAttribute = rv.getAttribute(Header.SAMPLE_NAME.getTitle());
                }
            }

            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public String getReseqTemplateID() {
            return reseqTemplateID;
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public void setPath(String path) {
            this.path = path;
        }
        
        public String getRGdata() {
            StringBuilder sb = new StringBuilder();
            String[] rgs = {this.rgId, this.rgPl, this.rgPu, this.rgLb, this.rgSm};
            for (int i = 0; i < rgs.length; i++) {
                if (i > 0 )
                    sb.append(":");
                sb.append(rgs[i]);
            }
            
            return sb.toString();
        }
    }
}
