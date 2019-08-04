#! /usr/local/bin/gawk -f
BEGIN {
	OFS = "\t"
	print "Sequence Name      |  H1-left | H1-right |  H2-left | H2-right | Read-Depth | Dup-Read |   in 1M | Col-read |   in 1M"
}

NR == FNR {
	cov[$1, $2] = $3
}

NR != FNR {
	meancov = 0.25 * (cov[$1, $2] + cov[$1, $3] + cov[$1, $4] + cov[$1, $5])
	if (meancov == 0) {
		#                                     | Read-Depth | Dup-Read |   in 1M | Col-read |   in 1M
		printf("%-18s | %8d | %8d | %8d | %8d |          - |        - |       - |        - |       -\n", $1, $2, $3, $4, $5)
	}
	else if ($6 == 0 && $7 == 0) {
		#                                            | Dup-Read |   in 1M | Col-read |   in 1M
		printf("%-18s | %8d | %8d | %8d | %8d | %10d |        - |       - |        - |       -\n", $1, $2, $3, $4, $5, meancov)
	}
	else if ($6 == 0) {
		#                                            | Dup-Read |   in 1M |
		printf("%-18s | %8d | %8d | %8d | %8d | %10d |        - |       - | %8d | %7.2f\n", $1, $2, $3, $4, $5, meancov, $7, $7 / meancov * 1000000)
	}
	else if ($7 == 0) {
		#                                                          | Col-Read |   in 1M
		printf("%-18s | %8d | %8d | %8d | %8d | %10d | %8d | %7.2f |        - |       -\n", $1, $2, $3, $4, $5, meancov, $6, $6 / meancov * 1000000)

	}
	else {
		printf("%-18s | %8d | %8d | %8d | %8d | %10d | %8d | %7.2f | %8d | %7.2f\n", $1, $2, $3, $4, $5, meancov, $6, $6 / meancov * 1000000, $7, $7 / meancov * 1000000)
	}
}