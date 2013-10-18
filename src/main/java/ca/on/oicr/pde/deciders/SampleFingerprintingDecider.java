/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;


import java.util.*;
import java.io.File;
import java.io.FileNotFoundException;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.Element;

/**
 *
 * @author pruzanov@oicr.on.ca
 */


public class SampleFingerprintingDecider extends OicrDecider {
    private Map<String, ReturnValue> pathToAttributes = new HashMap<String, ReturnValue>();
    
    private String output_prefix = "./";
    private String output_dir = "seqware-results";
    private String studyName;
    private String watchersList;
    
    
    //these params should come from settings xml file
    private String genomeFile;
    private String checkedSNPs;
    private int    checkPoints;
    
    //GATK params
    private double standCallConf = 50.0;
    private double standEmitConf = 10.0;
    private int    dcov          = 50;
    
    //Previous workflow runs
    private String existingMatrix;
    private String genotypes;
    
    private String currentRType; // resequencing type
    private String currentTType; // template type
    private String SNPConfigFile = "/.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml";

    public SampleFingerprintingDecider() {
        super();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("study-name", "Required: name of the study that we need to analyze.").withRequiredArg();
        parser.accepts("template-type", "Required: name of the study that we need to analyze.").withRequiredArg();
        parser.accepts("resequencing-type", "Optional: resequencing type for templates other than WG").withRequiredArg();
        parser.accepts("existing-matrix", "Optional: existing matrix from previous workflow run(s)").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to " 
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
	        + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("config-file","Optional. Path to a config file in .xml format "
                + "Default: /.mounts/labs/PDE/data/SampleFingerprinting/hotspots.config.xml").withRequiredArg();
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
        this.setMetaType(Arrays.asList("application/bam"));
	//Group by template type if no other grouping selected
        if (!this.options.has("group-by")) {
            this.setGroupingStrategy(Header.STUDY_SWA);
        } else {
            Log.warn("Passing group-by parameter overrides the default (STUDY_SWA) I hope you know what you are doing");
        }
        
        if (this.options.has("output-path")) {
          this.output_prefix = options.valueOf("output-path").toString();
          if (this.output_prefix.isEmpty())
	      this.output_prefix = "./";
          else if (!this.output_prefix.endsWith("/"))
              this.output_prefix += "/";
              
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
        
        if (this.options.has("watchers-list")) {
          String commaSepWatchers = options.valueOf("watchers-list").toString();
          
          String[] watchers = commaSepWatchers.split(",");
          for (String email : watchers) {
              if (email.contains("@oicr.on.ca")) {
                  this.watchersList += this.watchersList.isEmpty() ? email : "," + email;
              }
          }
          if (this.watchersList.contains("@"))
              this.watchersList = "";
	}
       
        if (this.options.has("study-name")) {
            this.studyName = options.valueOf("study-name").toString();
	} else {
            Log.error("study-name parameter needs to be set, aborting");
            System.exit(1);
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
            }
	}
        
        if (this.options.has("stand_emit_conf")) {
            try {
              this.standEmitConf = Double.valueOf(options.valueOf("stand_emit_call").toString());
            } catch (NumberFormatException nf) {
              Log.error("stand_emit_call failed to pass as double, make sure you supply a valid number");  
            }
	}
        
        if (this.options.has("dcov")) {
            try {
              this.dcov = Integer.valueOf(options.valueOf("dcov").toString());
            } catch (NumberFormatException nf) {
              Log.error("dcov failed to pass as int, make sure you supply a valid number");  
            }
	}
                               
        if (this.options.has("template-type")) {
          this.currentTType = options.valueOf("template-type").toString();
          this.currentRType = this.options.has("rsconfig-file") ? this.options.valueOf("rsconfig-file").toString() : "NA";
         
          boolean snpOK = this.configFromParsedXML(this.SNPConfigFile, this.currentTType, this.currentRType);
          if (!snpOK) {
              Log.error("Genotyping parameters not set up properly, aborting");
              System.exit(1);
          }
        } else {
              Log.error("template-type parameter needs to be set, aborting");
              System.exit(1);
        }
        
        //allows anything defined on the command line to override the defaults here.
        ReturnValue val = super.init();
        return val;
    }    


    protected String handleGroupByAttribute(String attribute, String group_id, String template_type) {

        String groupBy = this.getGroupingStrategy().getTitle();
        //Log.stdout("GROUP BY ATTRIBUTE: " + groupBy + " " + attribute);
        String[] parentNames = attribute.split(":");
        List <String> groupBySet = new ArrayList();
        groupBySet.add(parentNames[parentNames.length - 1]);
        if (null != template_type) {
            groupBySet.add(template_type);
        } 
               
        if (null != group_id) {
            groupBySet.add(group_id);
        }

        String[] myFilters = groupBySet.toArray(new String[0]);
        
        for (int i = 0; i < myFilters.length; i++) {
         if (null != myFilters[i]) {
          if (groupBy.length() > 1)
             groupBy = groupBy.concat(":" + myFilters[i]);
          else
             groupBy = groupBy.concat(myFilters[i]);
         }
        }

        return groupBy;
        
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
              
        Log.debug("CHECK FILE DETAILS:" + fm);
        //If there is no resequencing type, set it to NA
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String targetTemplateType     = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        if (null == targetResequencingType && !this.currentRType.equals("NA")) {
            return false;
        }
        if (null == targetTemplateType || !targetTemplateType.equals(this.currentTType)) {
            return false;
        }
              
	pathToAttributes.put(fm.getFilePath(), returnValue);
        return super.checkFileDetails(returnValue, fm);
    }
    
    
    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        //get files from study
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : vals) {
            String currVal  = r.getAttributes().get(groupBy);
            String group_id = r.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_group_id");
            String template_type = r.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            if (null == currVal || null == template_type || !template_type.equals(this.currentTType)) {
                continue;
            }
            currVal = null != group_id ? handleGroupByAttribute(currVal, group_id, template_type) : handleGroupByAttribute(currVal, null, template_type);
            
            if (null != currVal) {
                List<ReturnValue> vs = map.get(currVal);
                if (vs == null) {
                    vs = new ArrayList<ReturnValue>();
                }
                vs.add(r);
                map.put(currVal, vs);
            }
        }
     return map;
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);
        //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }
         
        //setting testmode for excluded stuff (I am not sure if it is necessery though, but just to be safe):
	if (commaSeparatedFilePaths.isEmpty()) {
	    this.setTest(true);
	}
        
        ReturnValue r = pathToAttributes.get(commaSeparatedFilePaths);
               
	Map<String, String> iniFileMap = new TreeMap<String, String>();
        iniFileMap.put("input_bams", commaSeparatedFilePaths);
        iniFileMap.put("output_prefix",this.output_prefix);
	iniFileMap.put("output_dir", this.output_dir);
            
        //TODO genotypes=
    
        if (null != this.checkedSNPs) {
          iniFileMap.put("checked_snps", this.checkedSNPs);
          iniFileMap.put("check_points", "" + this.checkPoints);
        }
        
        if (null != this.genomeFile) {
          iniFileMap.put("genome_file", this.genomeFile);  
        }
        
        if (null != this.genotypes && !this.genotypes.isEmpty()) {
          iniFileMap.put("genotypes", this.genotypes);
        } else {
          iniFileMap.put("genotypes", "");  
        }
        
        if (null != this.existingMatrix && !this.existingMatrix.isEmpty()) {
          iniFileMap.put("existing_matrix", this.existingMatrix);
        } else {
          iniFileMap.put("existing_matrix", "");  
        }

        if (null != this.watchersList && !this.watchersList.isEmpty()) {
          iniFileMap.put("watchers_list", this.watchersList);
        }
        
        return iniFileMap;
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
                                  this.checkedSNPs = cElement.getElementsByTagName("checked_snps").item(0).getTextContent();
                                  this.checkPoints = Integer.parseInt(cElement.getElementsByTagName("check_points").item(0).getTextContent());
                                  return true;
                              }                   
                    }   
                }	
	  }
        } catch(FileNotFoundException fnf) {
            System.err.println("File is not found");
        } catch (NullPointerException np) {
            Log.error("A value in config file for " + resequencingType + " is not initialized");
        } catch (NumberFormatException nf) {
            Log.error("A number of checkponts in the config file is not set properly, should be an integer");
        }  catch (Exception e) {
            e.printStackTrace();
        }
        return false;
      }
   
}
