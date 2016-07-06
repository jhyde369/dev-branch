#python separate_lock_train.py -s 1117329857 -e 1117429857 --min-seg-length 3600 --force --include-time-locked --include-time-until-end --threshold 25 --usertag "LowThresh" #small er7
#python separate_lock_train.py -s 1117329857 -e 1118336935 --min-seg-length 3600 --force --include-time-locked --include-time-until-end --threshold 25 --usertag "LowThresh" #full er7
#python separate_lock_train.py -s 1117329857 -e 1118336935 --min-seg-length 3600 --force --include-time-locked --include-time-until-end #full er7
#python plot_signif_trigs.py -s 1117329857 -e 1118336935 --min-seg-length 3600 # full er7

#python plot_signif_trigs.py -s 1117329857 -e 1117529857 --min-seg-length 3600 # full er7
#python separate_lock_train.py -s 1117329857 -e 1118336935 --min-seg-length 3600 --include-time-locked --include-time-until-end

#python split_data_train.py -s 1117329857 -e 1118336935 --min-seg-length 3600 --include-time-locked --include-time-until-end

python train_around_event.py -e 1126259462 --gracedb-ID "G184098" --include-time-locked --include-time-until-end #September event

python train_around_event.py -e 1128678900 --gracedb-ID "G197392" --include-time-locked --include-time-until-end #Other event

python train_around_event.py -e 1135136350 --gracedb-ID "G211117" --include-time-locked --include-time-until-end #Christmas event

