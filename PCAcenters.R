library(factoextra)
enriched <- read.csv(file.choose(), header = TRUE, row.names = 1)
enriched <- enriched[-c(48,50,56)]
colnames(enriched) <- c("AAA(K)","AAC(N)","AAG(K)","AAT(N)","ACA(T)","ACC(T)","ACG(T)","ACT(T)","AGA(R)","AGC(S)",
                        "AGG(R)","AGT(S)","ATA(I)","ATC(I)","ATT(I)","CAA(Q)","CAC(H)","CAG(Q)","CAT(H)","CCA(P)",
                        "CCC(P)","CCG(P)","CCT(P)","CGA(R)","CGC(R)","CGG(R)","CGT(R)","CTA(L)","CTC(L)","CTG(L)",
                        "CTT(L)","GAA(E)","GAC(D)","GAG(E)","GAT(D)","GCA(A)","GCC(A)","GCG(A)","GCT(A)","GGA(G)",
                        "GGC(G)","GGG(G)","GGT(G)","GTA(V)","GTC(V)","GTG(V)","GTT(V)","TAC(Y)","TAT(Y)","TCA(S)",
                        "TCC(S)","TCG(S)","TCT(S)","TGC(C)","TGT(C)","TTA(L)","TTC(F)","TTG(L)","TTT(F)")

enhanced <- read.csv(file.choose(), header = TRUE, row.names = 1)
enhanced <- enhanced[-c(48,50,56)]
colnames(enhanced) <- c("AAA(K)","AAC(N)","AAG(K)","AAT(N)","ACA(T)","ACC(T)","ACG(T)","ACT(T)","AGA(R)","AGC(S)",
                        "AGG(R)","AGT(S)","ATA(I)","ATC(I)","ATT(I)","CAA(Q)","CAC(H)","CAG(Q)","CAT(H)","CCA(P)",
                        "CCC(P)","CCG(P)","CCT(P)","CGA(R)","CGC(R)","CGG(R)","CGT(R)","CTA(L)","CTC(L)","CTG(L)",
                        "CTT(L)","GAA(E)","GAC(D)","GAG(E)","GAT(D)","GCA(A)","GCC(A)","GCG(A)","GCT(A)","GGA(G)",
                        "GGC(G)","GGG(G)","GGT(G)","GTA(V)","GTC(V)","GTG(V)","GTT(V)","TAC(Y)","TAT(Y)","TCA(S)",
                        "TCC(S)","TCG(S)","TCT(S)","TGC(C)","TGT(C)","TTA(L)","TTC(F)","TTG(L)","TTT(F)")

elevated <- read.csv(file.choose(), header = TRUE, row.names = 1)
elevated <- elevated[-c(48,50,56)]
colnames(elevated) <- c("AAA(K)","AAC(N)","AAG(K)","AAT(N)","ACA(T)","ACC(T)","ACG(T)","ACT(T)","AGA(R)","AGC(S)",
                        "AGG(R)","AGT(S)","ATA(I)","ATC(I)","ATT(I)","CAA(Q)","CAC(H)","CAG(Q)","CAT(H)","CCA(P)",
                        "CCC(P)","CCG(P)","CCT(P)","CGA(R)","CGC(R)","CGG(R)","CGT(R)","CTA(L)","CTC(L)","CTG(L)",
                        "CTT(L)","GAA(E)","GAC(D)","GAG(E)","GAT(D)","GCA(A)","GCC(A)","GCG(A)","GCT(A)","GGA(G)",
                        "GGC(G)","GGG(G)","GGT(G)","GTA(V)","GTC(V)","GTG(V)","GTT(V)","TAC(Y)","TAT(Y)","TCA(S)",
                        "TCC(S)","TCG(S)","TCT(S)","TGC(C)","TGT(C)","TTA(L)","TTC(F)","TTG(L)","TTT(F)")

enr.pca <- prcomp(enriched, scale = TRUE)
enh.pca <- prcomp(enhanced, scale = TRUE)
ele.pca <- prcomp(elevated, scale = TRUE)

fviz_pca_ind(enr.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping  
             )
fviz_pca_ind(enh.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping  
)
fviz_pca_ind(ele.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping  
)

enr.eigvec <- enr.pca$rotation[,1]
enh.eigvec <- enh.pca$rotation[,1]
ele.eigvec <- ele.pca$rotation[,1]

enr.eigvec
enh.eigvec
ele.eigvec
