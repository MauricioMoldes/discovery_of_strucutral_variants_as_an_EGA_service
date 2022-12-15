################
## Generate Sample Markdown
################

sed  "s/sample_template/$1/g" < sample_markdown_template.rmd > $1.rmd

