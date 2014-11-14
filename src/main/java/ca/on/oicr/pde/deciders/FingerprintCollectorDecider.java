/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
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

    private String outputPrefix = "./";
    private String outputDir = "seqware-results";
    private String manualOutput = "false";
    //GATK:
    private String gatkPrefix;
    private double standCallConf;
    private double standEmitConf;
    private int dcov;

    private String queue = "";
    private final static String BAM_METATYPE = "application/bam";

    private String templateTypeFilter = "";
    private String reseqTypeFilter = "";
    private String studyName;
    private Map<String, Map> reseqType;
    private final String SNPConfigFile = "/.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml";

    public FingerprintCollectorDecider () {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Optional. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
                + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (for example, to production)").withRequiredArg();
        parser.accepts("manual-output", "Optional. Set the manual output either to true or false "
                + "when running the workflow, the default is false").withRequiredArg();
        parser.accepts("gatk-memory", "Optional. Set the memory allocated to GATK jobs "
                + "when running the workflow, the default is 5000").withRequiredArg();
        parser.accepts("stand-call-conf", "Optional. Set GATK parameter stand_call_conf "
                + "when running the workflow, the default is 50.0").withRequiredArg();
        parser.accepts("stand-emit-conf", "Optional. Set GATK parameter stand_emit_conf "
                + "when running the workflow, the default is 10.0").withRequiredArg();
        parser.accepts("dcov", "Optional. Set GATK parameter dcov when running the workflow, the default is 200").withRequiredArg();
        parser.accepts("gatk-prefix", "Optional: the path to a dir on a low-latency filesystem for writing " 
                + "GATK temporary data. May prevent possible failures of a workflow run. Default: ./").withRequiredArg();

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(BAM_METATYPE));
        this.setGroupingStrategy(Header.FILE_SWA);
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
                this.dcov = Integer.valueOf(iniFileMap.get("dcov"));
            } else {
                Log.error("The given INI file does not exist: " + file.getAbsolutePath());
                System.exit(1);
            }

        }

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            this.templateTypeFilter = options.valueOf("template-type").toString();
            Log.debug("Setting template type is not necessary, however if set the decider will run the workflow only on this type of data");
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
        
        if (this.options.has("stand_call_conf")) {
            try {
                this.standCallConf = Double.valueOf(options.valueOf("stand_conf_call").toString());
            } catch (NumberFormatException nf) {
                Log.error("stand_conf_call failed to pass as double, make sure you supply a valid number");
                System.exit(1);
            }
        }

        if (this.options.has("stand_emit_conf")) {
            try {
                this.standEmitConf = Double.valueOf(options.valueOf("stand_emit_call").toString());
            } catch (NumberFormatException nf) {
                Log.error("stand_emit_call failed to pass as double, make sure you supply a valid number");
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

        if (this.options.has("verbose")) {
            Log.setVerbose(true);
        }

        if (this.options.has("manual-output")) {
            this.manualOutput = Boolean.valueOf(options.valueOf("manual-output").toString()) ? "true" : "false";
        }

        if (this.options.has("output-path")) {
            this.outputPrefix = options.valueOf("output-path").toString();
            if (!this.outputPrefix.endsWith("/")) {
                this.outputPrefix += "/";
            }
        }

        if (this.options.has("output-folder")) {
            this.outputDir = options.valueOf("output-folder").toString();
        }
        
        if (this.options.has("study-name")) {
            this.studyName = options.valueOf("study-name").toString();
	} else {
            Log.warn("study-name parameter is not set, will try to determine it automatically");
        }

        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String targetTemplateType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // If nulls, set to NA
        if (null == targetResequencingType || targetResequencingType.isEmpty()) {
            targetResequencingType = "NA";
        }
        if (null == targetTemplateType || targetTemplateType.isEmpty()) {
            targetTemplateType = "NA";
        }

        // Check filters
        if (!this.reseqTypeFilter.isEmpty() && !this.reseqTypeFilter.equals(targetResequencingType)) {
            return false;
        }
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equals(targetTemplateType)) {
            return false;
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
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    public ReturnValue customizeRun(WorkflowRun run) {

         //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }
        
        String inputFiles = "";
        String checkedSNPs = "";
        String checkPoints = "";
        
        for (FileAttributes atts : run.getFiles()) {
          if (checkedSNPs.isEmpty()) {
           String fileKey = fileSwaToSmall.get(atts.getOtherAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
           if (this.reseqType.containsKey(fileKey)) {
             checkedSNPs = this.reseqType.get(fileKey).get("file").toString();
             checkPoints = this.reseqType.get(fileKey).get("points").toString();
           } else {
                
           }
          }
          
          if (atts.getMetatype().equals(BAM_METATYPE)) {
            if (!inputFiles.isEmpty()) {
                inputFiles += ",";
            }
              inputFiles += atts.getPath();
          }
        }
                
        //setting testmode for excluded stuff (I am not sure if it is necessery though, but just to be safe):
	if (inputFiles.isEmpty()) {
	    this.setTest(true);
	}
        
        run.addProperty("input_files", inputFiles);
        run.addProperty("output_prefix",this.outputPrefix);
        run.addProperty("output_dir", this.outputDir);
        run.addProperty("gatk_prefix",this.gatkPrefix);
        run.addProperty("manual_output", this.manualOutput);
	
            
        if (null != checkedSNPs) {
          run.addProperty("checked_snps", checkedSNPs);
          run.addProperty("check_points", checkPoints);
        }
               
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

        
        return new ReturnValue();
    }

    private boolean configFromParsedXML(String fileName, String templateType, String resequencingType) {
        // SHOULD BE OK, though fields that need to be defined also need to be customized
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
        params.add(FingerprintCollectorDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
            if (null != fa.getLimsValue(Lims.GROUP_ID)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.GROUP_ID));
            }

            path = rv.getFiles().get(0).getFilePath() + "";
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

        public String getTissueType() {
            return tissueType;
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
    }
}
