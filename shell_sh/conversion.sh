for FILE in *.eps; do epspdf -b $FILE;temp=`expr "$FILE" | cut -d'.' -f1`; echo $temp.png;convert -quality 00 -density 100 -depth 3 $temp.pdf $temp.png;done
