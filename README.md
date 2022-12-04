To create tfs for specific run:  

	alien_find /alice/data/year/month/<your_run>/ *.tf | perl -p -e 's/^\/alice/alien:\/\/\/alice/' | tee tfs.lst




To fetch run  

	o2-raw-tf-reader-workflow --input-data tf.lst --onlyDet HMP -b  5| o2-hmpid-raw-to-digits-stream-workflow --fast-decode --ignore-dist-stf --get-results-statistics -b| o2-hmpid-digits-to-clusters-workflow --write-to-file -b --out-file hmpclus15.root



