# Each run uses one installed package version in a fresh R process.
# Usage: Rscript worker.R LIBRARY BACKEND OUTPUT.rds original|deterministic
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c(args[1], .libPaths()))
Sys.setenv(OMP_NUM_THREADS='1'); options(rf.cores=1)
suppressPackageStartupMessages(library(precmed))
data(countExample);data(survivalExample)
tau0 <- with(survivalExample,min(quantile(y[trt=='drug1'],.95),quantile(y[trt=='drug0'],.95)))
backend <- args[2]
profile <- if(length(args)>=4L)args[4] else 'original'
ns <- asNamespace('precmed')
original.boost <- if(backend=='gbm') getExportedValue('gbm','gbm') else getFromNamespace('fit_boosting','precmed')
record.limit <- as.integer(Sys.getenv('PRECMED_AUDIT_RECORD_LIMIT', '2'))
records <- list(); n.calls <- 0L
boost <- function(...) {
 a <- list(...)
 if(is.null(a$formula)) names(a)[1L] <- "formula"
 if(is.null(a$n.cores)) a$n.cores <- 1L
 before <- .Random.seed
 if(profile == 'deterministic') {
  caller.seed <- before
  on.exit(assign('.Random.seed',caller.seed,envir=globalenv()),add=TRUE)
  set.seed(100000L+n.calls)
  a$bag.fraction <- 1
  a$n.minobsinnode <- 50
  before <- .Random.seed
 }
 ans <- do.call(original.boost,a)
 n.calls <<- n.calls+1L
 if(length(records)<record.limit) {
  a$formula <- as.formula(paste(deparse(a$formula),collapse=' '),env=globalenv())
  iters <- unique(c(1L,10L,a$n.trees))
  pred <- suppressWarnings(predict(ans,a$data,n.trees=iters))
  best <- if(backend=='gbm') gbm::gbm.perf(ans,method='cv',plot.it=FALSE) else precmed:::best_boosting_iteration(ans,method='cv',plot.it=FALSE)
  records[[length(records)+1L]] <<- list(args=a,seed=before,predictions=pred,iterations=iters,best=best)
 }
 ans
}
for(name in c('intxcount','intxsurv')) {
 f <- get(name,ns)
 env <- new.env(parent=environment(f))
 assign(if(backend=='gbm') 'gbm' else 'fit_boosting',boost,env)
 environment(f) <- env
 unlockBinding(name,ns);assign(name,f,ns);lockBinding(name,ns)
}
specs <- readRDS(file.path(dirname(args[3]), 'specs.rds'))
results <- list()
saveRDS(results,args[3])
pdf(file=file.path(dirname(args[3]), paste0(backend,'-',profile,'-plots.pdf')))
for(i in seq_along(specs)) {
 spec <- specs[[i]]
 a <- lapply(as.list(spec$call)[-1],eval,envir=globalenv())
 response <- a$response; a$response <- NULL
 fun <- if(as.character(spec$call[[1]])=='catefit') {
   if(response=='count') catefitcount else catefitsurv
 } else if(response=='count') catecvcount else catecvsurv
 unsupported <- if(as.character(spec$call[[1]])=='catecv')setdiff(names(a),names(formals(catecv)))else character()
 # Survival CV has no dots: enforce one GBM core through the hook instead.
 if(response=='survival') a$n.cores <- NULL
 a$verbose <- 0
 a$plot.gbmperf <- FALSE
 cat('START',backend,spec$id,format(Sys.time()),'\n');flush.console()
 records <- list();n.calls <- 0L
 warnings <- character()
 elapsed <- system.time(result <- tryCatch(withCallingHandlers(do.call(fun,a),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity))
 cat('DONE',backend,spec$id,'elapsed',elapsed[['elapsed']],'boosters',n.calls,
     if(inherits(result,'error')) conditionMessage(result) else 'OK','\n');flush.console()
 results[[spec$id]] <- list(spec=spec,unsupported=unsupported,records=records,boost.calls=n.calls,
   elapsed=elapsed,warnings=unique(warnings),result=result,profile=profile,
   session=sessionInfo(),imports=packageDescription('precmed')$Imports)
 saveRDS(results,args[3])
}
dev.off()
