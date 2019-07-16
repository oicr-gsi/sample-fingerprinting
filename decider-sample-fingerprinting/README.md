## sample-fingerprinting Decider

Version 2.2.3

### Overview

This decider allows to customize and launch SampleFingerprinting workflow runs with the goal
to identify possible sample swaps.

### Compile

```
mvn clean install
```

### Usage
After compilation, [test](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#testing-the-workflow), [bundle](http://seqware.github.io/docs/3-getting-started/developer-tutorial/#packaging-the-workflow-into-a-workflow-bundle) and [install](http://seqware.github.io/docs/3-getting-started/admin-tutorial/#how-to-install-a-workflow) the workflow using the techniques described in the SeqWare documentation.

#### Options
These parameters can be overridden either in the INI file on on the command line using `--override` when [directly scheduling workflow runs](http://seqware.github.io/docs/3-getting-started/user-tutorial/#listing-available-workflows-and-their-parameters) (not using a decider). Defaults are in [square brackets].

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

    existing-matrix           string      Recycling results from previous SampleFingerprinting workflow runs
    config-file               string      specify alternative file with hotspot settings 
                                          (should not be needed)
    separate-platforms        true|false  Optional: Separate sequencer run platforms, determines if i.e.
                                          HiSeq and MiSeq data will be processed by separate runs 
                                          [true]
    after-date                string      Optional: Format YYYY-MM-DD. Implemented by OICRDecider, Only run
                                          on files modified after a certain date, not inclusive
    before-date               string      Optional: Format YYYY-MM-DD. Implemented by OICRDecider, Only run
                                          on files modified before a certain date, not inclusive
    watchers-list             string      List of people notified by email about detected sample swaps
    mixed-coverage            string      Parameter tells SampleFingerprinting workflow to
                                          use .fin files for calculation of similarities
                                          between the analyzed samples (disabled by default)
    queue                     string      Name of the (SGE) queue to schedule to 
                                          [empty]


### Support
For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .
