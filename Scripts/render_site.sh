#!/bin/bash
nohup Rscript -e "rmarkdown::render_site()"  > render_site.log &

