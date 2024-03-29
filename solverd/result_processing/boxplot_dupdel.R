b2 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch2_dupdel.txt", header=TRUE, sep="\t", dec=".")
b3 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch3_dupdel.txt", header=TRUE, sep="\t", dec=".")
b4 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch4_dupdel.txt", header=TRUE, sep="\t", dec=".")
b10 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch10_dupdel.txt", header=TRUE, sep="\t", dec=".")
b11 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch11_dupdel.txt", header=TRUE, sep="\t", dec=".")
b12 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch12_dupdel.txt", header=TRUE, sep="\t", dec=".")
b16 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch16_dupdel.txt", header=TRUE, sep="\t", dec=".")
b17 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch17_dupdel.txt", header=TRUE, sep="\t", dec=".")
b18 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch18_dupdel.txt", header=TRUE, sep="\t", dec=".")
b19 <- read.table("H:\\Data\\tmp3\\solverd\\documentation\\results\\dupdel_all\\batch19_dupdel.txt", header=TRUE, sep="\t", dec=".")


par(mfrow=c(1,2))
b2_dup <- boxplot(b2$Duplication, main="[ALL] Batch2 duplications", col="chartreuse4")
b2_del <- boxplot(b2$Deletion, main="[ALL] Batch2 deletions", col="firebrick4")
min(b2_dup$out)
max(b2_dup$out)
min(b2_del$out)
max(b2_del$out)

b3_dup <- boxplot(b3$Duplication, main="[ALL] Batch3 duplications", col="chartreuse4")
b3_del <- boxplot(b3$Deletion, main="[ALL] Batch3 deletions", col="firebrick4")
min(b3_dup$out)
max(b3_dup$out)
min(b3_del$out)
max(b3_del$out)

b4_dup <- boxplot(b4$Duplication, main="[ALL] Batch4 duplications", col="chartreuse4")
b4_del <- boxplot(b4$Deletion, main="[ALL] Batch4 deletions", col="firebrick4")
min(b4_dup$out)
max(b4_dup$out)
min(b4_del$out)
max(b4_del$out)

b10_dup <- boxplot(b10$Duplication, main="[ALL] Batch10 duplications", col="chartreuse4")
b10_del <- boxplot(b10$Deletion, main="[ALL] Batch10 deletions", col="firebrick4")
min(b10_dup$out)
max(b10_dup$out)
min(b10_del$out)
max(b10_del$out)

b11_dup <- boxplot(b11$Duplication, main="[ALL] Batch11 duplications", col="chartreuse4")
b11_del <- boxplot(b11$Deletion, main="[ALL] Batch11 deletions", col="firebrick4")
min(b11_dup$out)
max(b11_dup$out)
min(b11_del$out)
max(b11_del$out)

b12_dup <- boxplot(b12$Duplication, main="[ALL] Batch12 duplications", col="chartreuse4")
b12_del <- boxplot(b12$Deletion, main="[ALL] Batch12 deletions", col="firebrick4")
min(b12_dup$out)
max(b12_dup$out)
min(b12_del$out)
max(b12_del$out)

b16_dup <- boxplot(b16$Duplication, main="[ALL] Batch16 duplications", col="chartreuse4")
b16_del <- boxplot(b16$Deletion, main="[ALL] Batch16 deletions", col="firebrick4")
min(b16_dup$out)
max(b16_dup$out)
min(b16_del$out)
max(b16_del$out)

b17_dup <- boxplot(b17$Duplication, main="[ALL] Batch17 duplications", col="chartreuse4")
b17_del <- boxplot(b17$Deletion, main="[ALL] Batch17 deletions", col="firebrick4")
min(b17_dup$out)
max(b17_dup$out)
min(b17_del$out)
max(b17_del$out)

b18_dup <- boxplot(b18$Duplication, main="[ALL] Batch18 duplications", col="chartreuse4")
b18_del <- boxplot(b18$Deletion, main="[ALL] Batch18 deletions", col="firebrick4")
min(b18_dup$out)
max(b18_dup$out)
min(b18_del$out)
max(b18_del$out)

b19_dup <- boxplot(b19$Duplication, main="[ALL] Batch19 duplications", col="chartreuse4")
b19_del <- boxplot(b19$Deletion, main="[ALL] Batch19 deletions", col="firebrick4")
min(b19_dup$out)
max(b19_dup$out)
min(b19_del$out)
max(b19_del$out)

IQR()
