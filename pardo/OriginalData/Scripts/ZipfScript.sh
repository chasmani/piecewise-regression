#!/bin/sh
#If the console does not work type (for Mac): xattr -d com.apple.quarantine <YOUR FILE HERE>

####
#The idea of this script is to reformat Google's data to obtain a single table with all the words frequencies. This particular case is intended for german, but by changing every instance of ger for Google's corresponding abreviation of another language (e.g. eng for english) the table for any of the analysed languages can be obtained.
####

 
##This line will create one table per letter, eliminating words that have unwanted punctuation marks or numbers or numbers. For words whose frequency in one year is less than 50, it renders that particular frequency as 0 (This is done in order to eliminate pathological cases). 

for let in a b c d e f g h i j k l m n o p q r s t u v w x y z ; do awk '(!($1~"[0-9\!_\.]")){if($4>50) print $0}'  googlebooks-ger-all-1gram-20120701-$let > ALLger_$let; done

##Another way to do the same using awk itself.
#awk '(!($1~"[0-9\'\!_\.A-Z]")){if ($4>20) print $0}' googlebooks-eng-all-1gram-20120701-a #> ALL_aEng

##The first two brackets are regular expresions to eliminate some words, then a table is made, where the first column corresponds to the words, the next columns to the frequency of each word per year, and the last column to the total frequency of each word. 
printf '%s\0' ALL_* |xargs -0 cat | awk '{wordconc[$1]=$3; a[$1,$2]=$3;total[$1]+=$3} END{printf "word,"; for(y=1800;y<2012;++y){printf y",";} printf"total\r\n"; for (w in wordconc){ printf w","; for(y=1800;y<2010;++y) {if(a[w,y]=="") printf "0,"; else printf a[w,y]",";} printf total[w]"\r\n"} }' > Tables/all_words_per_year_ger.txt


