# Usage: Rscript inventory.R VIGNETTE_SOURCE_DIRECTORY OUTPUT_DIRECTORY
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
paths <- list.files(args[1],pattern='\\.Rmd$',full.names=TRUE)
specs <- inventory <- list()
for(path in paths) {
 lines <- readLines(path);inside <- FALSE; code<-character()
 for(line in lines) {
  if(grepl('^```\\{r',line)){header<-line;inside<-TRUE;code<-character()}
  else if(inside && line=='```') {
   inside<-FALSE
   if(!any(grepl('catefit\\(|catecv\\(',code))) next
   exprs <- parse(text=code)
   index<-0L
   for(expr in exprs) {
    call<-expr
    if(is.call(call) && as.character(call[[1]]) %in% c('<-','='))call<-call[[3]]
    if(!is.call(call) || !as.character(call[[1]]) %in% c('catefit','catecv'))next
    index<-index+1L
    uses <- grepl('boosting',paste(deparse(call),collapse=' '),fixed=TRUE)
    id<-paste(basename(path),sub('^```\\{r ([^,}]+).*','\\1',header),index,sep=':')
    info<-list(id=id,file=basename(path),header=header,uses.boosting=uses,call=call)
    inventory[[length(inventory)+1L]]<-info
    if(uses)specs[[length(specs)+1L]]<-info
   }
  } else if(inside)code<-c(code,line)
 }
}
saveRDS(specs,file.path(args[2], 'specs.rds'))
write.csv(do.call(rbind,lapply(inventory,function(s)data.frame(id=s$id,uses.boosting=s$uses.boosting,header=s$header))),file.path(args[2], 'inventory.csv'),row.names=FALSE)
cat('Identified',length(specs),'boosting scenarios\n')
