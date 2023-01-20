
# Read local CTF files:    
  o2-ctf-reader-workflow --ctf-input root.lst --max-tf 10 --delay 30   --onlyDet HMP -b | o2-hmpid-raw-to-digits-workflow -b --out-file digCTF_loc_10.root



o2-ctf-reader-workflow --ctf-input root.lst --max-tf 10 --delay 30   --onlyDet HMP -b | o2-hmpid-raw-to-digits--stream-workflow -b


# this writes file??
o2-ctf-reader-workflow --ctf-input root.lst --max-tf 10 --delay 30   --onlyDet HMP -b | o2-hmpid-raw-to-digits-workflow -b --out-file t.root | o2-hmpid-digits-to-root-workflow -b --in-file t.root --out-file hmpDig10_local_.root

