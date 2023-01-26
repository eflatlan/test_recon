
# Read local CTF files:    
  o2-ctf-reader-workflow --ctf-input root.lst --max-tf 128 --delay 20   --onlyDet HMP -b | o2-hmpid-raw-to-digits-workflow -b --out-file t.root |  o2-hmpid-digits-to-root-workflow -b --in-file t.root --out-file dig128.root



o2-ctf-reader-workflow --ctf-input root.lst --max-tf 10 --delay 30   --onlyDet HMP -b | o2-hmpid-raw-to-digits--stream-workflow -b


# this writes file??
o2-ctf-reader-workflow --ctf-input tfs.lst --max-tf 10 --delay 20   --onlyDet HMP -b | o2-hmpid-raw-to-digits-workflow -b --out-file t.root | o2-hmpid-digits-to-root-workflow -b --in-file t.root --out-file hmpDig10_remote_.root


o2-raw-tf-reader-workflow --input-data tfs_517623.lst --max-tf 400 --delay 20  --onlyDet HMP -b |
	o2-hmpid-raw-to-digits-stream-workflow --fast-decode --ignore-dist-stf --get-results-statistics -b |
 	o2-hmpid-digits-to-root-workflow -b --out-file tfs_517623_400tf_20delay.root

