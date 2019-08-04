#! /usr/local/bin/gawk -f

NR == FNR { # the first file is the signature summary file

	# copy the table to memory
	sign[FNR]["nm"] = $1
	sign[FNR][1] = $2 ; sign[FNR][2] = $3
	sign[FNR][3] = $4 ; sign[FNR][4] = $5
	sign[FNR]["it"] = $6 ; sign[FNR]["dt"] = $7
	# a reverse-search dictionary for signature's indeces
	pos[$1, 1, $2][FNR] ; pos[$1, 2, $3][FNR] ; pos[$1, 3, $4][FNR] ; pos[$1, 4, $5][FNR]

}

	# position indicator note:
	#           1  2           3  4
	#           |  |           |  |
	#    D site v  v I site    v  v
	# CCTCAGCCAGccgtGTTATAACTTAccgtTTACCAACTACATTTTTTGTAACGAACCAAA
	#           ^ I left clip  |  ^ I right clip
	#              |           ^ D left clip
	#              ^ D right clip

NR != FNR {

	if ($6 ~ /[SID]/) {
		split($6, seg, "[MSHID]", cig)
		# if left clipped
		if ($6 ~ /^[0-9]+S/) {
			if (length(pos[$3, 1, $4]) > 0) for (i in pos[$3, 1, $4]) if (index($10, sign[i]["it"])) {
				# print "read #" FNR, "signature #" i, "I LC"
				sign[i]["icount"] ++
			}
			if (length(pos[$3, 3, $4]) > 0) for (i in pos[$3, 3, $4]) if (index($10, sign[i]["dt"])) {
				# print "read #" FNR, "signature #" i, "D LC"
				sign[i]["dcount"] ++
			}
		}
		# if right clipped
		if ($6 ~ /[0-9]+S$/) {
			ep = $4 - 1
			for (j = 1 ; j <= length(cig) ; j++) if (cig[j] == "M" || cig[j] == "D") ep += seg[i]
			if (length(pos[$3, 2, ep]) > 0) for (i in pos[$3, 2, ep]) if (index($10, sign[i]["dt"])) {
				# print "read #" FNR, "signature #" i, "D RC"
				sign[i]["dcount"] ++
			}
			if (length(pos[$3, 4, ep]) > 0) for (i in pos[$3, 4, ep]) if (index($10, sign[i]["it"])) {
				# print "read #" FNR, "signature #" i, "I RC"
				sign[i]["icount"] ++
			}
		}
		# if insertion inside
		if ($6 ~ /[0-9]+I/) {
			ep = $4
			for (j = 1 ; j <= length(cig) ; j++) {
				if (cig[j] == "I") {
					if (length(pos[$3, 1, ep]) > 0) for (i in pos[$3, 1, ep]) { 
						if (sign[i][3] - sign[i][1] == seg[j] && index($10, sign[i]["it"])) {
							# print "read #" FNR, "signature #" i, "I Inside"
							sign[i]["icount"] ++
						}
					}
				}
				if (cig[j] == "M" || cig[j] == "D") ep += seg[j]
			}
		}
		# if deletion inside
		if ($6 ~ /[0-9]+D/) {
			ep = $4
			for (j = 1 ; j <= length(cig) ; j++) {
				if (cig[j] == "D") {
					if (length(pos[$3, 1, ep]) > 0) for (i in pos[$3, 1, ep]) {
						if (sign[i][3] - sign[i][1] == seg[j] && index($10, sign[i]["dt"])) {
							# print "read #" FNR, "signature #" i, "D Inside"
							sign[i]["dcount"] ++
						}
					}
				}
				if (cig[j] == "M" || cig[j] == "D") ep += seg[j]
			}
		}
	}
}

END {
	OFS = "\t"
	#OFMT = "%8.4f"
	for (i = 1 ; i <= length(sign) ; i++) {
		nm = sign[i]["nm"]
		print nm, sign[i][1], sign[i][2], sign[i][3], sign[i][4], sign[i]["icount"]+0, sign[i]["dcount"]+0
	}
}