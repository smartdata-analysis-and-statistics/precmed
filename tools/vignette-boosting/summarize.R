# Compare numeric scores, coefficients, validation curves, and ABC statistics.
# Usage: Rscript summarize.R GBM_RESULTS.rds GBM3_RESULTS.rds OUTPUT.csv
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
a<-readRDS(args[1])
b<-readRDS(args[2])
rows<-list()
for(id in intersect(names(a),names(b))) {
 x<-a[[id]]$result;y<-b[[id]]$result
 if(inherits(x,'error')||inherits(y,'error')) {
  rows[[length(rows)+1]]<-data.frame(scenario=id,component='ERROR',rmse=NA,maxabs=NA,exact=FALSE,error=if(inherits(x,'error'))conditionMessage(x)else conditionMessage(y));next
 }
 for(key in intersect(names(x),names(y))) {
  if(!grepl('^score\\.|^ate\\.|^coefficients$|^overall\\.ate',key))next
  u<-unlist(x[[key]]);v<-unlist(y[[key]])
  if(!is.numeric(u)||!is.numeric(v))next
  if(length(u)!=length(v))stop('Output shape mismatch: ',id,' ',key)
  ok<-is.finite(u)&is.finite(v)
  same.na<-identical(is.na(u),is.na(v)) && identical(is.infinite(u),is.infinite(v)) &&
    identical(sign(u[is.infinite(u)]),sign(v[is.infinite(v)]))
  if(!any(ok)) next
  rows[[length(rows)+1]]<-data.frame(scenario=id,component=key,rmse=sqrt(mean((u[ok]-v[ok])^2)),maxabs=max(abs(u[ok]-v[ok])),exact=same.na&&max(abs(u[ok]-v[ok]))<1e-8,error='')
 }
}
z<-do.call(rbind,rows);write.csv(z,args[3],row.names=FALSE)
print(z,row.names=FALSE)
