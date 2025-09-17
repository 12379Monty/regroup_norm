$css.css#!/bin/bash

css=CSS_FILE_NAME
group_normalization_documentation_md=group_normalization_documentation


readme_md=_README
 # readme
Rscript -e "rmarkdown::render('$readme_md.md', output_format = rmarkdown::html_document(),
    output_file='$readme_md.html')"  > $readme_md.log


# to generate htmls:


 # group_normalization_documentation
Rscript -e "rmarkdown::render('$group_normalization_documentation_md.md', output_format = rmarkdown::html_document(),
     output_file='$group_normalization_documentation_md.html')"  > $group_normalization_documentation_md.log









