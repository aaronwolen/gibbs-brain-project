# -t specify output format
# -s --standalone Produce output with an appropriate header and footer
# -S --smart produce typographically correct output
#   e.g., proper curly quotes, em dashes, etc.
# -i --incremental Make list items in slide shows display incremental

INPUT=lab-pres_methylation.md
OUTPUT=lab-pres_methylation_slides.html

pandoc -s -S -i -t dzslides --mathjax $INPUT -o $OUTPUT
