#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) stop("Need a input MSA file")

fa = readLines(args[1])
# maxmium number of sequences to keep
m_seq = 0
if (length(args) > 1) m_seq = as.numeric(args[2])

N_SEQ = length(fa)
N_POS = nchar(fa[2])

if (N_SEQ == 0) stop("Empty MSA file: ", args[1])

names = sub('.', '', sapply(strsplit(fa[seq(1, N_SEQ, 2)],"\\s+"), "[[", 1))
bases = strsplit(fa[seq(2, N_SEQ, 2)],"")

N_SEQ = length(names)
N_POS = nchar(fa[2])

seqs = matrix(1, ncol = N_POS, nrow = N_SEQ)

for(i in 1:N_SEQ) seqs[i, bases[[i]] == "-"] = 0

cetr = apply(seqs,2,sum) / dim(seqs)[1]

edist <- function(a, b = a, N = length(a)) { sum((a - b) ^ 2) / N; }

dists = apply(seqs, 1, edist, b = cetr, N = N_POS)

qs = quantile(dists, seq(0, 1, 0.25))

iqr = qs[4] - qs[2]

othres = qs[4] + 1.5 * iqr

cat("F_MSA: ", args[1], "\n", file = stderr())
cat("N_SEQ: ", N_SEQ, "\n", file = stderr())
cat("N_POS: ", N_POS, "\n", file = stderr())
cat("Q_000: ", qs[1], "\n", file = stderr())
cat("Q_025: ", qs[2], "\n", file = stderr())
cat("Q_050: ", qs[3], "\n", file = stderr())
cat("Q_075: ", qs[4], "\n", file = stderr())
cat("Q_100: ", qs[5], "\n", file = stderr())
cat(" Q_RM: ", othres, "\n", file = stderr())
cat(" N_RM: ", sum(dists > othres), "\n", file = stderr())

ord = order(dists)

dists = dists[ord]
names = names[ord]

n_seq = sum(dists <= othres)

if (m_seq > 0 && m_seq < n_seq) n_seq = m_seq

if (n_seq > 0) cat(names[1:n_seq], sep = '\n')

