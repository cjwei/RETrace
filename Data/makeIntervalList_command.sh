for f in ./*.info.txt
do
	fbname=$(basename "$f" .info.txt)
	./makeIntervalList.py --probes $f --prefix $fbname.temp
	cat IntervalList_header.txt $fbname.temp.targets.interval_list >$fbname.targets.interval_list
	cat IntervalList_header.txt $fbname.temp.probes.interval_list >$fbname.probes.interval_list
	rm $fbname.temp*
done
