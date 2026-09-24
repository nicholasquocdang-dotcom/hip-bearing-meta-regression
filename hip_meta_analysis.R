# Hip bearing meta-regression: main model, sensitivity analyses, and figures
# Requires: install.packages(c("metafor","clubSandwich")); hip_meta_data.csv in the working directory.
suppressMessages({library(metafor); library(clubSandwich)})
d <- read.csv("hip_meta_data.csv", stringsAsFactors=FALSE)
cc <- 5  # follow-up centered at 5 years
prep <- function(d, hr="hr", lo="lo", hi="hi", tcol="t"){
  d$yi <- log(d[[hr]]); d$sei <- (log(d[[hi]])-log(d[[lo]]))/(2*qnorm(.975)); d$vi <- d$sei^2
  # time term = (t - c) x [I(compared bearing is not MoXLPE) - I(reference bearing is not MoXLPE)]
  cmp <- sub(" vs .*","",d$contrast); ref <- sub(".* vs ","",d$contrast)
  d$time <- (d[[tcol]]-cc) * ((cmp!="MoXLPE") - (ref!="MoXLPE")); d }
d <- prep(d)
fit <- function(dd, mods){ f <- as.formula(paste("~ 0 +", paste(mods, collapse="+")))
  rma(yi, vi, mods=f, data=dd, method="REML", test="knha") }
tab <- function(r) round(cbind(HR=exp(coef(r)), lb=exp(r$ci.lb), ub=exp(r$ci.ub), p=r$pval),3)
het <- function(r) cat("k=",r$k," tau2=",round(r$tau2,4)," I2=",round(r$I2,1)," QE=",round(r$QE,2)," df=",r$k-r$p," QEp=",signif(r$QEp,3),"\n")
allm <- c("MoCPE","CoCPE","CoXLPE","CoC","MoM","time")
res <- fit(d, allm); cat("== MAIN MODEL (Hartung-Knapp) ==\n"); het(res); print(tab(res))
rob <- robust(res, cluster=d$registry, clubSandwich=TRUE)
cat("\n== MAIN MODEL (CR2 cluster-robust) ==\n"); print(round(cbind(HR=exp(coef(rob)), lb=exp(rob$ci.lb), ub=exp(rob$ci.ub), p=rob$pval, df=rob$dfs),3))
cat("\n== SENSITIVITY ==\n")
for (r in unique(d$registry)) { dd <- d[d$registry!=r,]; m <- allm[sapply(allm, function(v) any(dd[[v]]!=0))]
  cat("\nWithout", r, ": "); rr <- fit(dd, m); het(rr); print(tab(rr)) }
cat("\nExcluding MoM: "); rr <- fit(d[d$MoM==0,], setdiff(allm,"MoM")); het(rr); print(tab(rr))
cat("\n10-year estimates only: "); rr <- fit(d[d$t>=10,], setdiff(allm,"time")); het(rr); print(tab(rr))
d2 <- d; nj <- d2$registry=="NJR"; d2$hr[nj] <- d2$hr2[nj]; d2$lo[nj] <- d2$lo2[nj]; d2$hi[nj] <- d2$hi2[nj]; d2$t[nj] <- 2
d2 <- prep(d2); cat("\nNJR at 2 years instead of 10: "); rr <- fit(d2, allm); het(rr); print(tab(rr))
cat("\nAdjusted estimates only (drops EPRD): "); rr <- fit(d[d$adjusted=="yes",], allm); het(rr); print(tab(rr))
saveRDS(list(d=d,res=res,rob=rob), "fit.rds")

# ---- Figures ----
o <- readRDS("fit.rds"); d <- o$d; res <- o$res; b <- coef(res); vb <- vcov(res)
d$cmp <- sub(" vs .*","",d$contrast); d$ref <- sub(".* vs ","",d$contrast)
full <- c(MoCPE="Metal-on-conventional PE", CoCPE="Ceramic-on-conventional PE", CoXLPE="Ceramic-on-highly-cross-linked PE", CoC="Ceramic-on-ceramic", MoM="Metal-on-metal")
d$lab <- paste0(d$source, " (", d$registry, "): ", d$label)
ord <- names(full); vsx <- d$ref=="MoXLPE"; oth <- !vsx
tcrit <- qt(.975, df=res$k-res$p)
tot <- sum(vsx) + 3*length(ord) + sum(oth) + 3
rows<-c(); y<-c(); lb<-c(); ub<-c(); labs<-c(); hp<-c(); hl<-c(); pooled<-list(); cur <- tot
for (g in ord) { s <- d[vsx & d$cmp==g,]; hp<-c(hp,cur); hl<-c(hl,full[g]); cur<-cur-1
  for (j in seq_len(nrow(s))) { rows<-c(rows,cur); y<-c(y,s$yi[j]); lb<-c(lb,log(s$lo[j])); ub<-c(ub,log(s$hi[j])); labs<-c(labs,s$lab[j]); cur<-cur-1 }
  pooled[[g]] <- c(cur, b[g], sqrt(vb[g,g])); cur <- cur-2 }
hp<-c(hp,cur); hl<-c(hl,"Estimates against another reference bearing (HR as reported)"); cur<-cur-1
s <- d[oth,]; for (j in seq_len(nrow(s))) { rows<-c(rows,cur); y<-c(y,s$yi[j]); lb<-c(lb,log(s$lo[j])); ub<-c(ub,log(s$hi[j])); labs<-c(labs,s$lab[j]); cur<-cur-1 }
png("forest.png", width=2400, height=2900, res=230); par(mar=c(4,1,1,1))
fp <- forest(x=y, ci.lb=lb, ci.ub=ub, slab=labs, rows=rows, ylim=c(cur, tot+3), atransf=exp, at=log(c(0.25,0.5,1,2,4,8)),
  xlab="Hazard ratio for revision (log scale); groups are relative to metal-on-highly-cross-linked PE", header=c("Estimate","HR [95% CI]"), cex=0.6, refline=0, psize=0.9)
for (g in ord) { p <- pooled[[g]]; addpoly(x=p[2], sei=p[3], rows=p[1], mlab="Pooled (model, at 5 years)", atransf=exp, cex=0.6, efac=0.8, ci.lb=p[2]-tcrit*p[3], ci.ub=p[2]+tcrit*p[3]) }
text(fp$xlim[1], hp, hl, pos=4, font=2, cex=0.64); dev.off()
cols <- c(MoCPE="#1b6ca8", CoCPE="#f4a261", CoXLPE="#2a9d8f", CoC="#e76f51", MoM="#6d597a")
w <- 1/(d$vi + res$tau2)
png("bubble.png", width=2000, height=1400, res=220); par(mar=c(4.5,4.5,1,1))
plot(NA, xlim=c(0,11), ylim=log(c(0.45,6)), xlab="Follow-up time at which estimate applies (years)", ylab="HR vs metal-on-highly-cross-linked PE (log scale)", yaxt="n")
axis(2, at=log(c(0.5,1,2,4)), labels=c("0.5","1","2","4"), las=1); abline(h=0, lty=2, col="grey50")
set.seed(1)
for (g in ord) { s <- vsx & d$cmp==g; jit <- (seq_len(sum(s))-mean(seq_len(sum(s))))*0.06
  points(d$t[s]+jit, d$yi[s], pch=21, bg=adjustcolor(cols[g],0.55), col=cols[g], cex=0.6+2.2*sqrt(w[s]/max(w)))
  tt <- seq(0.5,10.5,0.1); lines(tt, b[g] + b["time"]*(tt-5), col=cols[g], lwd=2) }
legend("topleft", legend=full[ord], col=cols[ord], pt.bg=adjustcolor(cols[ord],0.55), pch=21, lwd=2, bty="n", cex=0.8); dev.off()
cat("ok\n")
