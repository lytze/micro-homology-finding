con <- file("data/cdsintrons.fa", "rt")
res_d <- "results/"
min_inter <- 10
max_inter <- 100
min_micro <- 5
gc_con <- 0.4

len <- data.frame(list(id = character(), len = integer(), count = integer()), stringsAsFactors = F)
n <- 0
repeat {
    ch <- readChar(con, 1)
    if (length(ch) == 0) {
        ## output
        len[n, ] <- list(gid, b, nrow(stack))
        if (nrow(stack) != 0) {
            sel <- nchar(gsub("[AT]", "", stack$Sequence)) / nchar(stack$Sequence) >= gc_con
            if (any(sel)) {
                res <- stack[sel, ]
                res$LenInt <- res[, 3] - res[, 2] - nchar(res[, 1])
                write.csv(res, paste0(res_d, gid, ".csv"), row.names = F)
            }
        }
        break
    }
    if (ch != "\n") { ## ignore new-lines
        if (ch == ">") {
            if (n > 0) {
                ## output
                len[n, ] <- list(gid, b, nrow(stack))
                if (nrow(stack) != 0) {
                    sel <- nchar(gsub("[AT]", "", stack$Sequence)) / nchar(stack$Sequence) >= gc_con
                    if (any(sel)) write.csv(stack[sel, ], paste0(res_d, gid, ".csv"), row.names = F)
                }
            }
            ## for each sequence
            n <- n + 1
            newchar <- ""
            gid <- ""
            while (newchar != "\n") {
                gid <- paste0(gid, newchar)
                newchar <- readChar(con, 1)
            }
            message("#", n, " ", gid)
            ## initialize
            dict <- new.env()
            stack <- data.frame(
                list(Sequence = character(), Start = integer(), End = integer()),
                stringsAsFactors = F)
            stack_top <- 1
            ## new sequence -> read the first N-mer
            b <- min_micro
            s <- readChar(con, min_micro)
        } else {
            ## not a new sequence -> move to the next N-mer
            b <- b + 1
            s <- paste0(substring(s, 2), ch)
        }
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
                            stack[stack_top - 1, ] <- list(paste0(stack$Sequence[stack_top - 1], ch),
                                                           dict[[s]], b)
                        } else {
                            ## merge and keep
                            stack[stack_top, ] <- list(paste0(stack$Sequence[stack_top - 1], ch),
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
    }
}
close(con)
write.csv(len, paste0(res_d, "stat.csv"), row.names = F)
