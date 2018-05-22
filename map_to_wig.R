require(parallel)
n_cores <- detectCores()

WI <- 50

res <- file("../results/mapped_miho_5.wig", "w")
res_d <- "../results/genes/"
cat("track type=wiggle_0 name=micro-homology-occurence-over-5\n", file = res)

N <- c("I", "II", "III")

for (ch in 1:3) {
    chr <- read.table(paste0("../data/chromosome", ch,".cds.coords"), stringsAsFactors = F)
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("chr", "res_d"))
    wig <- parLapplyLB(cl, 1:nrow(chr),
                       function(i) {
                           fi <- paste0(res_d, chr[i, 1], ".csv")
                           if (file.exists(fi)) {
                               miho <- read.csv(fi, header = T, stringsAsFactors = F)
                               map <- integer(max(miho[, 3]))
                               for (j in 1:nrow(miho)) {
                                   s <- miho[j, 2]:miho[j, 3]
                                   map[s] <- map[s] + 1
                               }
                               return(data.frame(loc = chr[i, 2]:(chr[i, 2] + length(map) - 1), map = map))
                           } else return(NULL)
                       })
    stopCluster(cl)
    wig_map <- integer(max(chr[, 3]))
    for (w in wig) {
        if (!is.null(w)) {
            wig_map[w$loc] <- ifelse(wig_map[w$loc] > w$map, wig_map[w$loc], w$map)
        }
    }
    cat(paste0("variableStep chrom=chr", N[ch], " span=", WI, "\n"), file = res, append = T)
    for (win in 0:(length(wig_map) %/% WI - 1)) {
        v <- max(wig_map[(win * WI + 1):((win + 1) * WI)])
        if (v != 0) cat(paste0(win * WI + 1, " ", v, "\n"), file = res, append = T)
    }
}

close(res)