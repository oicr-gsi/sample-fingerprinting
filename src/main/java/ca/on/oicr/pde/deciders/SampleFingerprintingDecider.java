/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;


import java.io.File;
import java.io.FileNotFoundException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
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
 *
 * @author pruzanov@oicr.on.ca
 */


public class SampleFingerprintingDecider extends OicrDecider {

    private SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private String output_prefix = "./";
    private String output_dir = "seqware-results";
    private String gatk_prefix = "./";
    private String studyName;
    private String watchersList = "";
    
    
    //these params should come from settings xml file
    private String genomeFile;
    
    //GATK params
    private double standCallConf = 50.0;
    private double standEmitConf = 10.0;
    private int    dcov          = 200;
    
    //Previous workflow runs and input files
    private String  existingMatrix;
    private String  genotypes;
    private String  inputFiles;
    private boolean provisionVcfs;
    

    private String templateTypeFilter = "";
    private String reseqTypeFilter = "";
    private Map<String, Map> reseqType;
    private String SNPConfigFile = "/.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml";
    private Map<String, BeSmall> fileSwaToSmall;
    
    public SampleFingerprintingDecider() {
        super();
        fileSwaToSmall  = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Optional: name of the study that we need to analyze.").withRequiredArg();
        parser.accepts("resequencing-type", "Optional: resequencing type for templates other than WG").withRequiredArg();
        parser.accepts("existing-matrix", "Optional: existing matrix from previous workflow run(s)").withRequiredArg();
        parser.accepts("provision-vcfs", "Optional: set to non-null re-using vcf files from the runs this decider will launch").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to " 
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
	        + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("config-file","Optional. Path to a config file in .xml format "
                + "Default: /.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml").withRequiredArg();
        parser.accepts("gatk-prefix", "Optional: the path to a dir on a low-latency filesystem for writing " 
                + "GATK temporary data. May prevent possible failures of a workflow run. Default: ./").withRequiredArg();
        parser.accepts("stand-call-conf", "Optional: GATK parameter, The minimum phred-scaled confidence threshold "
                + "at which variants should be called, Default: 50.0").withRequiredArg();
        parser.accepts("stand-emit-conf", "Optional: GATK parameter, Default: 10.0"
                + "The minimum phred-scaled confidence threshold at which variants should be emitted").withRequiredArg();
        parser.accepts("dcov", "Optional: GATK parameter, 50 for 4x, 200 for >30x WGS or Whole exome, Default: 50").withRequiredArg();
        parser.accepts("watchers-list", "Optional: Comma-separated list of oicr emails for people interested in monitoring this workflow").withRequiredArg();
    }
    
    @Override
    public ReturnValue init() {

        Log.debug("INIT");
        String [] metaTypes = {"application/bam","text/vcf-4","text/plain"}; 
        this.setMetaType(Arrays.asList(metaTypes));
        this.reseqType = new HashMap<String, Map>();
	//Group by template type if no other grouping selected
        if (!this.options.has("group-by")) {
            if (this.options.has("study-name")) {
                this.setGroupingStrategy(Header.STUDY_SWA);
            }
            if (this.options.has("root-sample-name")) {
                //TODO: net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header needs to be updated to support ROOT_SAMPLE_SWA
                //uncomment when fixed... this.setGroupingStrategy(Header.ROOT_SAMPLE_SWA);
                //error out until the above TODO is implemented
                throw new RuntimeException("ROOT_SAMPLE_SWA needs to be implemented in FindAllTheFiles.Header");
            }
            if (this.options.has("sequencer-run-name")) {
                this.setGroupingStrategy(Header.SEQUENCER_RUN_SWA);
            }
        } else {
            Log.warn("Passing group-by parameter overrides the defaults, I hope you know what you are doing");
        }
        
        if (this.options.has("output-path")) {
          this.output_prefix = options.valueOf("output-path").toString();
          if (this.output_prefix.isEmpty())
	      this.output_prefix = "./";
          else if (!this.output_prefix.endsWith("/"))
              this.output_prefix += "/";
              
        }
        
        if (this.options.has("gatk-prefix")) {
          this.gatk_prefix = options.valueOf("gatk-prefix").toString();
          if (this.gatk_prefix.isEmpty())
	      this.gatk_prefix = "./";
          else if (!this.gatk_prefix.endsWith("/"))
              this.gatk_prefix += "/";
        }
        
        if (this.options.has("output-folder")) {
          this.output_dir = options.valueOf("output-folder").toString();
          if (this.output_dir.isEmpty())
	      this.output_dir = "seqware-results";
	}
        
        if (this.options.has("existing-matrix")) {
          this.existingMatrix = options.valueOf("existing-matrix").toString();
          if (null == this.existingMatrix || this.existingMatrix.isEmpty())
	      this.existingMatrix = "";
	}
        
        if (this.options.has("provision-vcfs")) {
          this.provisionVcfs = options.valueOf("provision-vcfs").toString().isEmpty() ? false : true;
	}
        
        if (this.options.has("watchers-list")) {
          String commaSepWatchers = options.valueOf("watchers-list").toString();
          
          if (null != commaSepWatchers && !commaSepWatchers.isEmpty()) {
            Log.warn("We have " + commaSepWatchers + " for watchers");
            String[] watchers = commaSepWatchers.split(",");
            for (String email : watchers) {
              if (email.contains("@oicr.on.ca")) 
                  this.watchersList += this.watchersList.isEmpty() ? email : "," + email;
            }
          }
          if (this.watchersList.isEmpty() || !this.watchersList.contains("@"))
              this.watchersList = "";
	}
       
        if (this.options.has("study-name")) {
            this.studyName = options.valueOf("study-name").toString();
	} else {
            Log.warn("study-name parameter is not set, will try to determine it automatically");
        }
        
        if (this.options.has("config-file")) {
            this.SNPConfigFile = options.valueOf("config-file").toString();
            Log.warn("Got custom config file, will use settings from user's input");
	}
        
        if (this.options.has("genome-file")) {
            this.genomeFile = options.valueOf("genome-file").toString();
            Log.warn("Got genome file, will use it instead of the dafault one (hg19)");
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
                               
        if (this.options.has("template-type")) {
          this.templateTypeFilter = options.valueOf("template-type").toString();
        } else {
          Log.debug("template-type parameter ids not set, will include all available types from the study");
        }
        
        // Setting resequencing type works like setting up a filter, the decider should fiure out reseq. type automatically
        if (this.options.has("resequencing_type")) {
          this.reseqTypeFilter = this.options.valueOf("resequencing-type").toString();
        } else {
          Log.debug("resequencing_type parameter ids not set, will include all available types from the study");
        }
        
        
        if (this.options.has("rsconfig-file")) {
          this.SNPConfigFile = this.options.valueOf("rsconfig-file").toString();
        }       
        
        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }    
    
    

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
              
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String targetTemplateType     = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        // If nulls, set to NA
        if (null == targetResequencingType || targetResequencingType.isEmpty())
            targetResequencingType = "NA";
        if (null == targetTemplateType || targetTemplateType.isEmpty())
            targetTemplateType = "NA";

        // Check filters
        if (!this.reseqTypeFilter.isEmpty() && !this.reseqTypeFilter.equals(targetResequencingType))
            return false;
        if (!this.templateTypeFilter.isEmpty() && !this.templateTypeFilter.equals(targetTemplateType))
            return false;
        
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
            if (!myTypes.contains(r.getFiles().get(0).getMetaType().toString()))
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
    public ReturnValue customizeRun(WorkflowRun run) {
        Log.debug("INI FILE:" + run.getIniFile().toString());
        //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }
        
        Set<String> vcfFiles = new HashSet<String>();
        this.inputFiles = "";
        this.genotypes  = "";
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
          
          if (atts.getMetatype().equals("application/bam")) {
            if (!this.inputFiles.isEmpty()) {
                this.inputFiles += ",";
            }
              this.inputFiles += atts.getPath();
          } else if (atts.getMetatype().equals("text/vcf-4")) { // Add them to set, sort later
              vcfFiles.add(atts.getPath());
          } else if (atts.getMetatype().equals("text/plain")) { // There should be only one
              this.existingMatrix = atts.getPath();
          }
        }
        
        // Match vcf files to bam files
        boolean vcfFoundOnce = false;
        for (String bam : this.inputFiles.split(",")) {
            String bamBase = bam.substring(bam.lastIndexOf("/") + 1,bam.lastIndexOf(".bam"));
            boolean vcfFound = false;
            for (String vcfPath : vcfFiles) {
                if (vcfPath.contains(bamBase)) {
                   if (!this.genotypes.isEmpty()) {
                     this.genotypes += ",";
                   }
                   this.genotypes += vcfPath; 
                   vcfFound = true;
                   vcfFoundOnce = true;
                }
            }
            this.genotypes += vcfFound ? "" : ",";
        }
        
        if (!vcfFoundOnce) {
            this.genotypes = " ";
        }
         
        //setting testmode for excluded stuff (I am not sure if it is necessery though, but just to be safe):
	if (this.inputFiles.isEmpty()) {
	    this.setTest(true);
	}
        
        run.addProperty("input_files", this.inputFiles);
        run.addProperty("output_prefix",this.output_prefix);
        run.addProperty("gatk_prefix",this.gatk_prefix);
	run.addProperty("output_dir", this.output_dir);
            
        if (null != checkedSNPs) {
          run.addProperty("checked_snps", checkedSNPs);
          run.addProperty("check_points", checkPoints);
        }
        
        if (null != this.genomeFile) {
          run.addProperty("genome_file", this.genomeFile);  
        } 
        
        run.addProperty("stand_call_conf", "" + this.standCallConf);
        run.addProperty("stand_emit_conf", "" + this.standEmitConf);
        run.addProperty("dcov", "" + this.dcov);
        
        if (null != this.genotypes && !this.genotypes.isEmpty()) {
          run.addProperty("genotypes", this.genotypes);
        } // genotypes set to 'space' if there are no vcfs, so the case when this is supposed to be empty handled already
        
        run.addProperty("study_name", this.studyName);
                
        if (null != this.existingMatrix && !this.existingMatrix.isEmpty()) {
          run.addProperty("existing_matrix", this.existingMatrix);
        } else {
          run.addProperty("existing_matrix", " ");
        }
        
        if (this.provisionVcfs) {
          run.addProperty("provision_vcfs", "TRUE");  
        } else {
          run.addProperty("provision_vcfs", " ");
        }
        
        if (null != this.watchersList && !this.watchersList.isEmpty()) {
          run.addProperty("watchers_list", this.watchersList);
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
                    if (!templateType.equals(nElement.getAttribute("id")))
                        continue;
                }
                
                if (nNode.hasChildNodes()) {
                    NodeList children = nNode.getChildNodes();
                    for (int cld = 0; cld < children.getLength(); cld++) {
                        Node cNode = children.item(cld);
                              if (cNode.getNodeType() == Node.ELEMENT_NODE) {
                                  Element cElement = (Element) cNode;
                                  if (!resequencingType.equals(cElement.getAttribute("id")))
                                      continue;
                                  Map<String, String> rtypeData = new HashMap<String, String>();
                                  rtypeData.put("file", cElement.getElementsByTagName("checked_snps").item(0).getTextContent());
                                  rtypeData.put("points", cElement.getElementsByTagName("check_points").item(0).getTextContent());
                                  this.reseqType.put(templateType.concat(resequencingType), rtypeData);
                                  return true;
                              }                   
                    }   
                }	
	  }
        } catch(FileNotFoundException fnf) {
            Log.error("File is not found");
        } catch (NullPointerException np) {
            Log.error("A value in config file for " + resequencingType + " is not initialized");
        } catch (NumberFormatException nf) {
            Log.error("A number of checkponts in the config file is not set properly, should be an integer");
        }  catch (Exception e) {
            e.printStackTrace();
        }
        return false;
      }
    
    public static void main(String args[]){
 
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(SampleFingerprintingDecider.class.getCanonicalName());
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

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            // Having metatype as part of details is needed since we deal with multiple mime types her
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode() + fa.getMetatype();
            // We are going to group by template type only (if we did not receive template type as a parameter)
            groupByAttribute = fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);
            String tgtReseq  = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
            if (null != tgtReseq) {
                groupByAttribute = groupByAttribute.concat(tgtReseq);
            } else {
                groupByAttribute = groupByAttribute.concat("NA");
            }
            path = rv.getFiles().get(0).getFilePath().toString();
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
    }
   
}
