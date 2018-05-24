require(seqinr)
require(parallel)
n_cores <- detectCores()

fa <- read.fasta("../data/artif.fa")
res_d <- "../results/artif/"
min_inter <- 10
max_inter <- 100
min_micro <- 5
gc_con <- 0.4

cl <- makeCluster(n_cores)
clusterExport(cl, c("len", "res_d", "min_inter", "max_inter", "min_micro", "gc_con"))
len <- parSapplyLB(cl, fa, function(seq) {
    gid <- attr(seq, "name")
    ## initialize
    dict <- new.env()
    stack <- data.frame(
        list(Sequence = character(), Start = integer(), End = integer()),
        stringsAsFactors = F)
    stack_top <- 1
    ## new sequence -> read the first N-mer
    b <- min_micro
    s <- paste(seq[1:b], collapse = "")
    while (b <= length(seq)) {
        if (is.null(dict[[s]])) {
            ## N-mer not recorded yet -> record the position
            dict[[s]] <- b
        } else if (b - dict[[s]] >= min_inter) {
            ## spacer is larger than mininal value ->
            if (b - dict[[s]] <= max_inter) {
                ## spacer is good ->
                if (stack_top == 1 || stack$End[stack_top - 1] != b - 1) {
                    ## empty stack || no adjacent N-mer -> push to result stack
                    stack[stack_top, ] <- list(s, dict[[s]], b)
                    stack_top <- stack_top + 1
                } else {
                    ## adjacent N-mer on stack top -> check and merge
                    stt <- which(stack$Sequence == stack$Sequence[stack_top - 1] &
                                     stack$Start == dict[[s]] - 1)
                    if (length(stt)) {
                        ## need for merge ->
                        if (stt == stack_top - 1) {
                            ## merge and substitute
                            stack[stack_top - 1, ] <- list(paste0(stack$Sequence[stack_top - 1], seq[b - 1]),
                                                           dict[[s]], b)
                        } else {
                            ## merge and keep
                            stack[stack_top, ] <- list(paste0(stack$Sequence[stack_top - 1], seq[b - 1]),
                                                       stack$Start[stt] + 1, b)
                            stack_top <- stack_top + 1
                        }
                    } else {
                        stack[stack_top, ] <- list(s, dict[[s]], b)
                        stack_top <- stack_top + 1
                    }
                }
            }
            dict[[s]] <- b
            ## requirement for refresh dict: 1. empty, 2. spacer longer than min
        }
        ## move to next N-mer
        b <- b + 1
        s <- paste0(seq[(b - min_micro + 1):b], collapse = "")
    }
    ## output
    if (nrow(stack) != 0) {
        sel <- nchar(gsub("[at]", "", stack$Sequence)) / nchar(stack$Sequence) >= gc_con
        if (any(sel)) write.csv(stack[sel, ], paste0(res_d, gid, ".csv"), row.names = F)
    }
    return(c(gid, b, nrow(stack)))
})
stopCluster(cl)
writeLines(c("id,len,count", apply(len, 2, paste, collapse = ",")),
           paste0(res_d, "stat.csv"))

