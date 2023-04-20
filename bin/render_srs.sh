#!/usr/bin/env bash

Rscript -e "rmarkdown::render('srs_curve.rmd', output_file='$PWD/srs_curve.pdf', output_format='pdf_document', clean=TRUE, knit_root_dir='$PWD')"
Rscript -e "rmarkdown::render('srs_curve.rmd', output_file='$PWD/srs_curve.html', output_format='html_document', clean=TRUE, knit_root_dir='$PWD')"
