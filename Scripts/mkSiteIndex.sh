#!/bin/bash


# Makes two fimes: _site.yml and _index.Rmd 
# index.Rmd needs to be updated each time a new html file
# - typically a rendered Rmd script - is added to the folder


OUTPUT="_site.yml"
TITLE="Training and Data Re-Analysis"

DESCRIPTION="Section: ZZZ"

echo 'name: "reAnalysis"' > $OUTPUT
echo 'output_dir: "."' >> $OUTPUT
echo 'navbar:' >> $OUTPUT
echo '  title: '$TITLE >> $OUTPUT
echo '  left:' >> $OUTPUT
echo '    - text: "Home"' >> $OUTPUT
echo '      href: index.html' >> $OUTPUT

OUTPUT="_index.Rmd"    ####  will build index.Rmd

echo '---' > $OUTPUT
echo 'title: '$TITLE >> $OUTPUT
echo '---' >> $OUTPUT
echo >> $OUTPUT
echo $DESCRIPTION >> $OUTPUT
echo >> $OUTPUT

files="`ls *.html`"
for f in $files
do
  if [ $f != "index.html" ]
  then
    echo '<li><a href="'$f'">'${f}'</a></li>' >> $OUTPUT
  fi
done

###     echo '<li><a href="'$f'">'${f#_}'</a></li>' >> $OUTPUT

