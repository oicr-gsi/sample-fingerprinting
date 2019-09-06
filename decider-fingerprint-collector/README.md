## fingerprint-collector Decider

Version 1.3.3

### Overview

Fingerprint Collector Decider pre-configures workflow runs that accepts alignments ( .bam files )
from Novoalign (or BWA) aligners. The purpose is to provide a "lighter" version of analysis that would be
capable of processing large-scale data more efficiently than 1.0 version. Fingerprint Collector workflow run would
produce .vcf.gz, .tbi and .fin files for downstream processing (Sample Fingerprinting 2.0 workflow).


### Compile

```
mvn clean install
```

### Usage
After compilation, [test](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#testing-the-workflow), [bundle](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#packaging-the-workflow-into-a-workflow-bundle) and [install](http://seqware.github.io/docs/3-getting-started/admin-tutorial/#how-to-install-a-workflow) the workflow using the techniques described in the SeqWare documentation.

#### Options
These parameters can be overridden either in the INI file or on the command line using `--override` when [directly scheduling workflow runs](http://seqware.github.io/docs/3-getting-started/user-tutorial/#listing-available-workflows-and-their-parameters) (not using a decider). Defaults are in [square brackets].

Required:

    study-name                string      A required parameter passed by the decider
                                          or on the command line if workflow is launched
                                          manually
    template-type             string      type of experimental template - WG (Whole Genome sequencing)
                                          EX (exome) etc. (works as a filter)
    resequencing-type         string      resequencing type - also a filter, specifies which resequencing
                                          type to process

Input/output:

    output-prefix             dir         The root output directory
    output-dir                string      The sub-directory of output_prefix where 
                                          the output files will be moved
    manual-output             true|false  When false, a random integer will be 
                                          inserted into the path of the final file 
                                          in order to ensure uniqueness. When true,
                                          the output files will be moved to the 
                                          location of output_prefix/output_dir
                                          [false]

Optional:

    preprocess-bam            true|false  Flag that indicates if re-ordering/adding Read Groups
                                          are needed, should rarely be true
                                          [false]
    gatk-prefix               string      prefix for customizing temporary files' location for GATK
    gatk-memory               integer     RAM in Mb allocated to GATK jobs      
    queue                     string      Name of the (SGE) queue to schedule to [production]


### Support
For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .
