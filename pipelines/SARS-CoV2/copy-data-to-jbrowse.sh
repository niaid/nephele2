#!/bin/bash
pair_ids=$1
sample_ids=$2
fcid=$3
grouping_regex=$4

echo "pair_ids: $pair_ids"
echo "sample_ids: $sample_ids"
echo "fcid: $fcid"
echo "grouping_regex: $grouping_regex"

# Create a new trackList for this fcid
# update_trackList.py called below will look
# for this file and append tracks to it
echo '{"formatVersion": 1, "tracks": []}' > ${fcid}_trackList.json

# For pairs (ex: A + B, unmerged), just update the track info
# Data will be sent altogether into the sample_id dir (next)
for pair_id in ${pair_ids}
do
	# pair_id might have a comma, or leading
	# or trailing square bracket, get rid of it
	pair_id=$(echo $pair_id | sed "s/[][,]//g")

	# use the --bw flag to only create the bw track
	# send the grouping_regex to get the sample_id
	# from the pair_id (for the category)
	update_trackList.py \
		--fcid="${fcid}" \
		--id="${pair_id}" \
		--bw \
		--regex="${grouping_regex}"
done

for sample_id in ${sample_ids}
do
	# sample_id might have a comma, or leading
	# or trailing square bracket, get rid of it
	sample_id=$(echo $sample_id | sed "s/[][,]//g")

	# create the dir where we will send all the data for this
	# sample_id then rsync all the data for this sample_id there
	dir=/usr/share/nginx/html/jbrowse/data/tracks/samples/${fcid}/${sample_id}
	ssh jbrowse@covid-19.bio.nyu.edu "mkdir -p ${dir}"
	rsync --copy-links ${sample_id}[-_]* \
		jbrowse@covid-19.bio.nyu.edu:${dir}/.

	# update the trackList for all sample_id (merged) files (default)
	update_trackList.py --fcid="$fcid" --id="$sample_id"
done

