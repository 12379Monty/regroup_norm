#!/bin/bash
##

css=CSS_FILE_NAME
GN4GD_implementationGuide_md=GN4GD_implementationGuide
article_summary_md=article_summary
TerminologyVariabillityOmic_md=TerminologyVariabillityOmic


readme_md=_README
 # readme
Rscript -e "rmarkdown::render('$readme_md.md', output_format = rmarkdown::html_document(),
    output_file='$readme_md.html')"  > $readme_md.log

# to generate htmls:
 ################################

 # GN4GD_analysisGuide
Rscript -e "rmarkdown::render('$GN4GD_implementationGuide_md.md', output_format = rmarkdown::html_document(),
    output_file='$GN4GD_implementationGuide_md.html')"  > $GN4GD_implementationGuide_md.log
     
 
# article_summary
Rscript -e "rmarkdown::render('$article_summary_md.md', output_format = rmarkdown::html_document(),
     output_file='$article_summary_md.html')"  > $article_summary_md.log


# TerminologyVariabillityOmic
Rscript -e "rmarkdown::render('$TerminologyVariabillityOmic_md.md', output_format = rmarkdown::html_document(),
     output_file='$TerminologyVariabillityOmic_md.html')"  > $TerminologyVariabillityOmic_md.log

