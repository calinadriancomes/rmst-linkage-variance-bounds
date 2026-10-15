# Optional source-side demonstration. Requires R and the official survival 3.8-12 archive.
# Supply the archive as argument 1 and a new output directory as argument 2.
main <- function() {
args <- commandArgs(TRUE)
if (length(args) != 2) stop('Usage: Rscript produce_marginal_tables.R survival_3.8-12.tar.gz output_dir')
folder <- tempfile('survival_'); dir.create(folder)
on.exit(unlink(folder, recursive=TRUE))
untar(args[1], exdir=folder, tar='internal')
pkg <- file.path(folder,'survival')
stopifnot(read.dcf(file.path(pkg,'DESCRIPTION'))[1,'Version'] == '3.8-12')
load(file.path(pkg,'data','cancer.rda'))
z <- rotterdam; tau <- 1826.25; n <- nrow(z)
stopifnot(n==2982, length(unique(z$pid))==n, !anyNA(z))
rec <- z$recur==1
in_window_death <- !rec & z$death==1 & z$dtime<=z$rtime
late <- !rec & z$death==1 & z$dtime>z$rtime
pt <- ifelse(rec,z$rtime,ifelse(in_window_death,z$dtime,z$rtime))
pe <- as.integer(rec | in_window_death)
dt <- z$dtime; de <- z$death
stopifnot(sum(rec)==1518,sum(in_window_death)==152,sum(late)==43,sum(pe==0)==1312)
stopifnot(sum(pt>=tau)==1571,sum(dt>=tau)==2084)
stopifnot(all(de==0 | ifelse(pe==1,pt<=dt,pt<dt)))
make_table <- function(t,e) {
  ev <- sort(unique(t[e==1 & t<tau]))
  data.frame(time=ev, n_risk=vapply(ev,function(s)sum(t>=s),integer(1)),
             n_event=vapply(ev,function(s)sum(t==s & e==1),integer(1)))
}
ptab <- make_table(pt,pe);dtab <- make_table(dt,de)
stopifnot(nrow(ptab)==878,nrow(dtab)==600)
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
write.csv(ptab,file.path(args[2],'P_event_risk.csv'),row.names=FALSE,quote=FALSE)
write.csv(dtab,file.path(args[2],'D_event_risk.csv'),row.names=FALSE,quote=FALSE)
message('Supply the accompanying metadata.json, which fixes definitions and checkpoints.')

}
main()
