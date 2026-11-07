# Second check of rigid_permutations() and patch_permutations(), with node_test() itself and a
# generator not used to build them: a field of exponential covariance, an L-shaped region of
# random sites, sparse grids and a narrow strip, 150 replicates, 199 permutations.
# Run from the root of the package, with the package installed. About 15 min on 10 cores.
library(comdynamite)
cores <- as.integer(Sys.getenv("CORES", "10"))
n_rep <- as.integer(Sys.getenv("N_REP","150")); n_perm <- 199
disc <- function(u, c0, r) sqrt(rowSums(sweep(u, 2, c0)^2)) < r
lshape <- function(m) { repeat { u <- cbind(runif(3*m,0,24), runif(3*m,0,24)); u <- u[u[,1] < 12 | u[,2] < 12, ]; if (nrow(u) >= m) return(u[1:m, ]) } }
grid <- as.matrix(expand.grid(1:20, 1:20))
land <- function(s) switch(s,
  off_adj = { u <- lshape(280); list(u, disc(u, c(5,6), 3.5)*1 - disc(u, c(5,15), 3.5)*1) },
  off_rand = { u <- lshape(280); repeat { c1 <- u[sample(280,1),]; c2 <- u[sample(280,1),]; a <- disc(u,c1,3); b <- disc(u,c2,3); if (!any(a&b) && sum(a)>=8 && sum(b)>=8) break }; list(u, a*1-b*1) },
  sp50_adj = { u <- grid[sort(sample(400,200)),]; a <- u[,1]%in%6:10 & u[,2]%in%8:12; b <- u[,1]%in%11:15 & u[,2]%in%8:12; list(u, a*1-b*1) },
  sp30_mid = { u <- grid[sort(sample(400,120)),]; a <- u[,1]%in%3:8 & u[,2]%in%8:13; b <- u[,1]%in%12:17 & u[,2]%in%8:13; list(u, a*1-b*1) },
  narrow = { u <- as.matrix(expand.grid(1:60, 1:3)); a <- u[,1]%in%10:17; b <- u[,1]%in%22:29; list(u, a*1-b*1) })
one <- function(i, s, r) { set.seed(7000 + i); L <- land(s); x <- unname(L[[1]]); ev <- L[[2]]; m <- nrow(x)
  if (min(sum(ev==1),sum(ev==-1)) < 5) return(NULL)
  rownames(x) <- paste0("s", seq_len(m))
  U <- chol(exp(-as.matrix(dist(x))/r) + diag(1e-8, m)); h <- drop(crossprod(U, rnorm(m))); h <- (h-min(h))/diff(range(h))*100
  tree <- treeshapesim::sim_tree_shape(12); tr <- species_traits(tree, lambda=1, optima_range=c(0,100), rmax_range=c(0,30))
  comm <- simulate_static(tr, h, incidence=TRUE)$abundance; comm <- comm[,colSums(comm)>0,drop=FALSE]; if (ncol(comm)<4) return(NULL)
  rownames(comm) <- rownames(x); nodes <- node_contrasts(ape::keep.tip(tree,colnames(comm)))
  b <- cbind(event=ev); rownames(b) <- rownames(x); env <- cbind(noise=rnorm(m)); rownames(env) <- rownames(x)
  f <- function(p) { if (is.null(p) || nrow(p$event) < 19) return(NA_real_); q <- node_test(comm,nodes,env,b,n_perm=n_perm,area_permutations=p)$biogeography$p_partial_site; mean(q[,"event"]<=.05,na.rm=TRUE) }
  pr <- tryCatch(suppressWarnings(rigid_permutations(x,b,n_perm)), error=function(e) NULL); pp <- patch_permutations(x,b,n_perm)
  c(rigid=f(pr), patch=f(pp), layouts=if (is.null(pr)) 0 else nrow(pr$event)) }
results <- list()
for (s in c("off_adj","off_rand","sp50_adj","sp30_mid","narrow")) for (r in c(4)) {
  o <- parallel::mclapply(seq_len(n_rep), function(i) tryCatch(one(i,s,r),error=function(e) conditionMessage(e)), mc.cores=cores, mc.preschedule=FALSE)
  bad <- vapply(o,is.character,TRUE); if (any(bad)) print(table(unlist(o[bad])))
  m <- do.call(rbind,o[!bad & !vapply(o,is.null,TRUE)]); st <- function(v) { v <- v[!is.na(v)]; sprintf("%.1f (%.1f) n=%d",100*mean(v),100*sd(v)/sqrt(length(v)),length(v)) }
  results[[s]] <- colMeans(m, na.rm = TRUE)
  cat(s,"rigid",st(m[,"rigid"]),"| patch",st(m[,"patch"]),"| layouts",round(mean(m[,"layouts"])),"\n") }
saveRDS(do.call(rbind, results), "data-raw/validation_nulls_check.rds")
