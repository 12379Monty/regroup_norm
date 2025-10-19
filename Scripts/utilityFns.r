 # Save and load
 ##################################################
 # save single object with name ObjName  to file FileName
 saveObj <- function(FileName='', ObjName='', DataDir=file.path(WRKDIR, 'Data')){
  assign(FileName, get(ObjName))
  save(list=FileName, file=file.path(DataDir, FileName))
  rm(list=FileName)
 }

 # load single object stored in FileName and assign to ObjName in local env
 loadObj <- function(FileName='', ObjName='', DataDir=file.path(WRKDIR, 'Data')){
  load(file.path(DataDir, FileName))
  assign(ObjName, get(FileName),pos=1)
  rm(list=FileName)
 }

 # Copy def to help_DIR
 #dput(saveObj, file.path(help_DIR, 'saveObj.r'))
 #dput(loadObj, file.path(help_DIR, 'loadObj.r'))

 # timing 
 ##################################################
 startTimedMessage <- function(...) {
        x <- paste0(..., collapse='')
        message(x, appendLF=FALSE)
        ptm <- proc.time()
        return(ptm)
 }
 stopTimedMessage <- function(ptm) {
        time <- proc.time() - ptm
        message(" ", round(time[3],2), "s")
 }

 # help
 ################################################
 static_help <- function(pkg, topic, out, links = tools::findHTMLlinks()) {
  pkgRdDB = tools:::fetchRdDB(file.path(find.package(pkg), 'help', pkg))
  force(links)
  tools::Rd2HTML(pkgRdDB[[topic]], out, package = pkg,
                 Links = links, no_links = is.null(links))
 }

 # kappa for clusering aggreement with labels
 ################################################
 getKappa <- function(Clustering.v, Labels.v)
  {
     require(concord)
     # Reorder if names are availableo
     if(!is.null(names(Clustering.v)))
     Labels.v <- Labels.v[names(Clustering.v)]

     MS_1.v  <- as.numeric(as.character(factor(Labels.v,
                  level=c('MSS', 'MSI'), labels=c(1,2))))
     MS_2.v  <- as.numeric(as.character(factor(Labels.v,
                  level=c('MSS', 'MSI'), labels=c(2,1))))
     kappa.1 <- cohen.kappa(cbind(MS_1.v, Clustering.v))
     kappa.2 <- cohen.kappa(cbind(MS_2.v, Clustering.v))
     max(kappa.1$kappa.c,kappa.2$kappa.c)

     kappa.v <- round(max(kappa.1$kappa.c,kappa.2$kappa.c),2)
  }
 #dput(getKappa, file.path(help_DIR, 'getKappa.r'))


 # legend
 twoRowLegend <- function(Legend.vec, Loc='top', Nrow=2, ...) {
  # draws horizontal legend in 1 or 2 rows

  if(Nrow==2) if(length(Legend.vec) %% 2) Legend.vec <- c(Legend.vec,NA)

  legend_order <- matrix(1:length(Legend.vec), nrow=Nrow, byrow=F)

  legend(Loc, horiz=F, lty=1, bty='n',
         legend=names(Legend.vec)[legend_order],
         col=Legend.vec[legend_order], lwd=2,
     ncol=ncol(legend_order),...)
}


 ##################################
 panel.hist <- function(x, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(usr[1:2], 0, 1.5) )
         h <- hist(x, plot = FALSE)
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
         legend('bottomright', paste0('IQR=',round(IQR(x),2)), bty='n', text.col='blue')
     }
##############################
kablePlus_f <- 
function (Table, caption = NULL, monospace = FALSE, striped = TRUE, 
    full_width = FALSE, format = c("html", "latex"), caption_color = "black") 
{
   # Table: to be diplayed (may have to be data.frame)
   # the rest are optional parameters  

    if (nrow(Table) == 0) {
        return(NULL)
    }
    format <- match.arg(format)
    if (striped) {
        bootstrap_options <- "striped"
    }
    else {
        bootstrap_options <- "basic"
    }
    caption <- kableExtra::text_spec(caption, format, color = caption_color, 
        bold = TRUE)
    kable <- Table %>% kableExtra::kable(caption = caption) %>% 
        kableExtra::kable_styling(bootstrap_options, full_width = full_width) %>% 
        kableExtra::row_spec(0:nrow(Table), extra_css = "border-bottom: 1px solid;") %>% 
        kableExtra::column_spec(1:ncol(Table), border_left = TRUE, 
            border_right = TRUE, monospace = monospace, color = "black")
    if (monospace) {
        kable <- kable %>% kableExtra::row_spec(0, monospace = TRUE)
    }
    return(kable)
}
