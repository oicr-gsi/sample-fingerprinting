##sample-fingerprinting workflow

Version 2.0.2

###Overview

Sample Fingerprinting workflow produces a graphical report on genotype clustering by sample (donor) and sends alerts if there's a possible swap/mix-up detected. It uses .vcf .tbi and .fin files produced by FingerprintCollector workflow. The below graph describes the process:

![sample-fingerprinting flowchart](docs/SampleFingerprinting_specs.png)

###Dependencies

This workflow requires:

* [SeqWare](http://seqware.github.io/)
* [tabix](http://sourceforge.net/projects/samtools/files/tabix/) 0.2.6
* [vcftools] (http://vcftools.sourceforge.net/) 0.1.10

###Compile

```
mvn clean install
```

###Usage
After compilation, [test](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#testing-the-workflow), [bundle](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#packaging-the-workflow-into-a-workflow-bundle) and [install](http://seqware.github.io/docs/3-getting-started/admin-tutorial/#how-to-install-a-workflow) the workflow using the techniques described in the SeqWare documentation.

####Options
These parameters can be overridden either in the INI file on on the command line using `--override` when [directly scheduling workflow runs](http://seqware.github.io/docs/3-getting-started/user-tutorial/#listing-available-workflows-and-their-parameters) (not using a decider). Defaults are in [square brackets].

Required:

    study_name                string      A required parameter passed by the decider
                                          or on the command line if workflow is launched
                                          manually

Input/output:

    output_prefix             dir         The root output directory
    output_dir                string      The sub-directory of output_prefix where 
                                          the output files will be moved
    manual_output             true|false  When false, a random integer will be 
                                          inserted into the path of the final file 
                                          in order to ensure uniqueness. When true,
                                          the output files will be moved to the 
                                          location of output_prefix/output_dir
                                          [false]

Optional:

    data_dir                  dir         A standard SeqWare parameter specifying the
                                          sub-directory where the output files will 
                                          be staged
    checked_snps              string      path to the vcf file (bgzipped and tabix-indexed)
                                          that contains coordinates and ids of dbSNPs
                                          used to genotype .bam files
    check_points              int         Number of snps in 'checked_snps' file
    exisiting_matrix          string      A file with previously calculated jaccard indexes 
                                          (recycled from previous runs to speed up calculations)
                                          Should be passed by the decider if available
    watchers_list             string      A comma-delimited list of emails 
                                          (of those who are interested in monitoring 
                                          SampleFingerprinting workflow runs)
    queue                     string      Name of the (SGE) queue to schedule to [production]



###Output files
**sample_fingerprint.[StudyName].report.zip**
Contains index.html with Sample Swap report, similarity matrix files, heatmaps of clustered samples in png format

###Support
For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .
