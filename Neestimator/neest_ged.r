library(dartR)

ged <- readRDS("./data/ged.rda")
table(pop(ged))
jerra <- gl.keep.pop(ged, pop.list="Jerrabomberra West")

jerra2 <- gl.filter.callrate(jerra, method="loc", threshold = 1)
jerra3 <- gl.filter.allna(jerra2)

res <- gl.LDNe(jerra3,neest.path = "d:/programms/Neestimator/", outpath = "d:/temp", outfile = "jw.txt")
