# Replay representative booster inputs, not entire workflows, under common settings.
# Usage: Rscript replay.R BASELINE_RESULTS.rds GBM3_LIBRARY OUTPUT.csv
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
.libPaths(c(args[2],.libPaths()))
suppressPackageStartupMessages(requireNamespace('gbm',quietly=TRUE));requireNamespace('gbm3',quietly=TRUE)
fit_boosting<-getFromNamespace('fit_boosting','precmed')
best_boosting_iteration<-getFromNamespace('best_boosting_iteration','precmed')
reference <- readRDS(args[1])
rows <- list()
for(case in reference) {
 if(inherits(case$result,'error') || !length(case$records)) next
 for(j in seq_along(case$records)) {
  record<-case$records[[j]]
  for(profile in c('vignette','fixed_trees','full_bag','full_bag_min50','full_bag_min100','adapter_vs_compat')) {
   a<-record$args
   if(is.null(a$shrinkage))a$shrinkage<-.1
   if(profile %in% c('fixed_trees','full_bag','full_bag_min50','full_bag_min100')) a$cv.folds<-0
   if(profile %in% c('full_bag','full_bag_min50','full_bag_min100'))a$bag.fraction<-1
   if(profile=='full_bag_min50')a$n.minobsinnode<-50
   if(profile=='full_bag_min100')a$n.minobsinnode<-100
   .Random.seed<-record$seed
   invisible(capture.output(old<-do.call(if(profile=='adapter_vs_compat') gbm3::gbm else gbm::gbm,a)))
   .Random.seed<-record$seed
   b<-a
   if(b$cv.folds==0)b$fold.id<-rep(1L,nrow(b$data))
   new<-do.call(fit_boosting,b)
   if(profile %in% c('vignette','adapter_vs_compat')) {
    old.best<-max(10,if(profile=='adapter_vs_compat') gbm3::gbm.perf(old,method='cv',plot.it=FALSE) else gbm::gbm.perf(old,method='cv',plot.it=FALSE))
    new.best<-max(10,best_boosting_iteration(new,method='cv',plot.it=FALSE))
    counts<-list(c(old.best,new.best))
   } else counts<-list(c(10,10),c(a$n.trees,a$n.trees))
   for(count in counts) {
    p<-suppressWarnings(predict(old,a$data,n.trees=count[1]))
    q<-suppressWarnings(predict(new,a$data,n.trees=count[2]))
    rows[[length(rows)+1]]<-data.frame(scenario=case$spec$id,fit=j,distribution=a$distribution,profile,old.trees=count[1],new.trees=count[2],correlation=cor(p,q),rmse=sqrt(mean((p-q)^2)),maxabs=max(abs(p-q)),exact=max(abs(p-q))<1e-8)
   }
  }
 }
}
z<-do.call(rbind,rows);write.csv(z,args[3],row.names=FALSE)
print(z,row.names=FALSE)
