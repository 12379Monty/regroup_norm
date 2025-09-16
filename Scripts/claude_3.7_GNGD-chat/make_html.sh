$css.css#!/bin/bash

css=CSS_FILE_NAME
affy_tiling_arrays_md=affy_tiling_arrays
citations_summary_md=citations_summary
group_normalization_summary_md=group_normalization_summary
TerminologyVariabillityOmic_md=TerminologyVariabillityOmic


readme_md=_README
 # readme
Rscript -e "rmarkdown::render('$readme_md.md', output_format = rmarkdown::html_document(),
    output_file='$readme_md.html')"  > $readme_md.log


# to generate htmls:

 # GN4GD_analysisGuide
Rscript -e "rmarkdown::render('$affy_tiling_arrays_md.md', output_format = rmarkdown::html_document(),
    output_file='$affy_tiling_arrays_md.html')"  > $affy_tiling_arrays_md.log
     

 # citations_summary
Rscript -e "rmarkdown::render('$citations_summary_md.md', output_format = rmarkdown::html_document(),
     output_file='$citations_summary_md.html')"  > $citations_summary_md.log


 # group_normalization_summary
Rscript -e "rmarkdown::render('$group_normalization_summary_md.md', output_format = rmarkdown::html_document(),
     output_file='$group_normalization_summary_md.html')"  > $group_normalization_summary_md.log



# TerminologyVariabillityOmic
Rscript -e "rmarkdown::render('$TerminologyVariabillityOmic_md.md', output_format = rmarkdown::html_document(),
     output_file='$TerminologyVariabillityOmic_md.html')"  > $TerminologyVariabillityOmic_md.log







