## Overview

R code used for the machine learning portion of: "Machine learning approaches for influenza A risk assessment: identifying predictive correlates using ferret model in vivo data". Manuscript under review at Communications Biology.

This project also includes CSV files for summary data & metrics used to create figures from the manuscript in the "inputs" directory. 

Study makes use of previously published data:
## References and Resources
Kieran TJ, Sun X, Creager HM, Tumpey TM, Maines TR, Belser JA. An aggregated dataset of serial morbidity and titer measurements from influenza A virus-infected ferrets. Accepted Sci Data.

Kieran TJ, Sun X, Maines TR, Beauchemin CAA, Belser JA. 2024. Exploring associations between viral titer measurements and disease outcomes in ferrets inoculated with 125 contemporary influenza A viruses. J Virol 98:e01661-23. (PMID 38240592)
([https://doi.org/10.1128/jvi.01661-23](https://doi.org/10.1128/jvi.01661-23))

## Note
R script uses predicted sialic acid binding preference (RBS) and predicted polymerase activity (PBS) markers (based on molecular sequence), along with selected HA and PB2 gene markers, that are not included in the CSV file from Sci Data, but may be cross-referenced as reported in Kieran et al (PMID 38240592). 

## Manuscript Abstract
In vivo assessments of influenza A virus (IAV) pathogenicity and transmissibility in ferrets represent a crucial component of many pandemic risk assessment rubrics, but few systematic efforts to identify which data from in vivo experimentation are most useful for predicting pathogenesis and transmission outcomes have been conducted. To this aim, we aggregated viral and molecular data from 125 contemporary IAV (H1, H2, H3, H5, H7, and H9 subtypes) evaluated in ferrets under a consistent protocol. Three overarching predictive classification outcomes (lethality, morbidity, transmissibility) were constructed using machine learning (ML) techniques, employing datasets emphasizing virological and clinical parameters from inoculated ferrets, limited to viral sequence-based information, or combining both data types. Among 11 different ML algorithms tested and assessed, gradient boosting machines and random forest algorithms yielded the highest performance, with models for lethality and transmission consistently better performing than models predicting morbidity. Comparisons of feature selection among models was performed, and highest performing models were validated with results from external risk assessment studies. Our findings show that ML algorithms can be used to summarize complex in vivo experimental work into succinct summaries that inform and enhance risk assessment criteria for pandemic preparedness that take in vivo data into account. 

##
##
##
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

