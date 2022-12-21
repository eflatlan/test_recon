To create tfs for specific run:  

	alien_find /alice/data/year/month/<your_run>/ *.tf | perl -p -e 's/^\/alice/alien:\/\/\/alice/' | tee tfs.lst



Load TFs from grid locally   

	alien_cp -T 1 /alice/data/2022/JUN/517623/raw/1250 file://

Where -T <number> specifies number of jobs to run in paralell
Will load subframe 1250 of TFs in run 517623 in a file named 1250 in working directory.




To fetch run  

	o2-raw-tf-reader-workflow --input-data tf.lst --onlyDet HMP -b  5 |  
	o2-hmpid-raw-to-digits-stream-workflow --fast-decode --ignore-dist-stf --get-results-statistics -b|  
	o2-hmpid-digits-to-clusters-workflow --write-to-file -b --out-file hmpclus15.root

Digit  


	o2-raw-tf-reader-workflow --input-data tfs_529005.lst --max-tf 180 --delay 200  --onlyDet HMP -b  5 |
	o2-hmpid-raw-to-digits-stream-workflow --fast-decode --ignore-dist-stf --get-results-statistics -b |
 	o2-hmpid-digits-to-root-workflow -b --out-file hmp529005_Dig_180.root



	o2-raw-tf-reader-workflow --input-data tfs_529691.lst --max-tf 120 --delay 800  --onlyDet HMP -b |
	o2-hmpid-raw-to-digits-stream-workflow --fast-decode --ignore-dist-stf --get-results-statistics -b |
 	o2-hmpid-digits-to-root-workflow -b --out-file hmp529691_Dig_120.root


LHC22t/529691

LHC22q/529005


alien_find /alice/data/year/LHC22t/529691/ *.tf | perl -p -e 's/^\/alice/alien:\/\/\/alice/' | tee tfs_529691.lst




