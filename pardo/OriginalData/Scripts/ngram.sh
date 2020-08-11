#! /bin/bash
#
# Author: Jay McKinnon <jay@opendna.com>
# URL: http://opendna.com/ngram
# License: CC-BY-NC 
# http://creativecommons.org/licenses/by-nc/3.0/
# Ref: http://books.google.com/ngrams/datasets
#
#
##########################################################
#  WARNING  WARNING  WARNING  WARNING  WARNING  WARNING  #
##########################################################
# It is really, REALLY easy to fill every last sector of #
# your hard drive with this script. Be cautious and pay  #
# attention.                                             #
#                                                        #
# This script will download ten ~210mb ZIP files, one at #
# a time, unpack the ~1000mb CSV inside, and (grep)      #
# search for each keyword. It will write the output to a #
# keyword-named CSV, and then delete the source file.    #
#                                                        #
# You must have a MINIMUM 1gb to run the default script. #
# You may remove the lines to delete source files, so    #
# can run more searches without downloading anew, but    #
# you must have upwards of 10gb available. If you change #
# the script to process 2-grams or higher, WATCH OUT!    #
# A multi-keyword search of 5-grams without deletes can  #
# easily top a terabyte of data!                         #
##########################################################

for count in 0 1 2 3 4 5 6 7 8 9
# starts a loop for these ngram files
  do wget http://commondatastorage.googleapis.com/books/ngrams/books/googlebooks-eng-all-1gram-20090715-$count.csv.zip
  # download each ngram file
    mv googlebooks-eng-all-1gram-20090715-$count.csv.zip* 1gram-$count.zip
    # rename downloaded ZIP file to something manageable
    unzip 1gram-$count.zip
    # unpack the ZIP archive
    mv googlebooks-eng-all-1gram-20090715-$count.csv 1gram-$count.csv
    # rename the unpacked CSV to something manageable
    for word in telegraph telephone television Internet internet
    # starts a subloop for your search terms
      do grep $word 1gram-$count.csv >> $word.csv
      # appends each keyword hit to the end of the keyword file
    done
    # ends subloop
     rm 1gram-$count.*
     # deletes ZIP and CSV files
     # Remove to save bandwith and fill your hard drive.
done
# ends mainloop
end
# kills script
##########################################################
