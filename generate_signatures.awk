#!/usr/local/bin/gawk -f
# usage
# <SZ=[feature size]> <TH=[feature diff threshold]> reference.fa mh.1.out <mh.2.out ...>

BEGIN {
	if (SZ == "") SZ = 10 ; if (TH == "") TH = 2
	printf("Feature size:          %7d\n", SZ) | "cat >& 2"
	printf("Specificity threshold: %7d\n", TH) | "cat >& 2"
	print "-------------------------------------------" | "cat >& 2"
	nname = 0;
	iname = 0;
	tot = 0 ; outrange = 0 ; fail = 0
}

NR == FNR {
	if ($0 ~ "^>") {
		nname++
		sname[nname] = substr($1, 2)
		block = 0
	}
	else {
		if (block == 0) g[sname[nname], "b"] = length($0)
		g[sname[nname], block] = $0
		block += length($0)
		g[sname[nname], "l"] = block
	}
}

function g_extract(name, s, e,
	ss, ee, i, b, r) {
	if (g[name, "l"] == "") return - 1
	if (s < 1) return -1
	if (e > g[name, "l"]) return -1
	if (s > e) return -1
	b = g[name, "b"]
	s--
	e--
	ss = s - (s % b)
	ee = e - (e % b)
	# print name, b, s, ss, e, ee
	if (ss == ee) {
		r = substr(g[name, ss], s % b + 1, e - s + 1)
	}
	else {
		r = substr(g[name, ss], s % b + 1)
		for (i = ss + b ; i < ee ; i += b) {
			r = r g[name, i]
		}
		r = r substr(g[name, ee], 1, e % b + 1)
	}
	return r
}

function gen_signature(n, s, e, k, r,
	k1l, k1r, k2l, k2r, mh, dis, i) {
	k1l = g_extract(n, s-SZ, s-1)
	k1r = g_extract(n, s+k, s+k+SZ)
	k2l = g_extract(n, e-SZ, e-1)
	k2r = g_extract(n, e+k, e+k+SZ)
	mh = g_extract(n, s, s+k-1)
	# print k1l, g_extract(n, s, s+k-1), k1r
	# print k2l, g_extract(n, e, e+k-1), k2r
	if (k1l == -1 || k1r == -1 || k2l == -1 || k2r == -1) return -1
	if (k1l == k2l || k1r == k2r) return -2
	dis = 0
	for (i = 1 ; i <= SZ ; i ++) {
		if (substr(k1l, i, 1) == substr(k2l, i, 1)) dis ++
	}
	if (dis < TH) return -2
	dis = 0
	for (i = 1 ; i <= SZ ; i ++) {
		if (substr(k1r, i, 1) == substr(k2r, i, 1)) dis ++
	}
	if (dis < TH) return -2
	r[1] = k2l mh k1r
	r[2] = k1l mh k2r
	return 0
}

NR != FNR {
	if (FNR == 1) iname ++
	# print SZ, TH, SN, $1, $2, $3
	tot ++
	test = gen_signature(sname[iname], $1, $2, $3, sign)
	if (test == 0) {
		print sname[iname], $1, $1+$3-1, $2, $2+$3-1, sign[1], sign[2]
		# the output records are:
		# sequence name | I/D clip sites 1 2 3 4 |
		# WT signature 1 | WT signature 2 | I signature | D signature
		# position indicator note:
		#           1  2           3  4
		#           |  |           |  |
		#    D site v  v I site    v  v
		# CCTCAGCCAGccgtGTTATAACTTAccgtTTACCAACTACATTTTTTGTAACGAACCAAA
		#           ^ I left clip  |  ^ I right clip
		#              |           ^ D left clip
		#              ^ D right clip
	}
	else if (test == -1) outrange ++
	else if (test == -2) fail ++
}

END {
	printf("Total MH pairs:        %7d\n", tot) | "cat >& 2"
	printf("Successed:             %7d (%.1f%)\n", tot - outrange - fail, (tot - outrange - fail) / tot * 100)| "cat >& 2"
	printf("Signature out of range:%7d (%.1f%)\n", outrange, outrange / tot * 100) | "cat >& 2"
	printf("Specificity failed:    %7d (%.1f%)\n", fail, fail / tot * 100) | "cat >& 2"
}