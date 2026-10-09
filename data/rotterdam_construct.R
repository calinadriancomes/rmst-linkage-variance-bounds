# Construct recorded endpoints from survival 3.8-12 without fitting any model.
if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args <- commandArgs(TRUE)
if (length(args) != 2) stop("Usage: Rscript data/rotterdam_construct.R survival_3.8-12.tar.gz rotterdam_check.csv")
folder <- tempfile(); dir.create(folder)
untar(args[1], exdir=folder, tar="internal")
desc <- read.dcf(file.path(folder,"survival","DESCRIPTION"))
stopifnot(desc[1,"Version"] == "3.8-12")
load(file.path(folder,"survival","data","cancer.rda"))
z <- rotterdam; n <- nrow(z); tau <- 1826.25
stopifnot(n==2982, length(unique(z$pid))==n, !anyNA(z), all(z$recur %in% 0:1), all(z$death %in% 0:1), all(z$rtime>=0), all(z$dtime>=0))
stopifnot(!any(z$recur==1 & z$death==1 & z$dtime<z$rtime))
rec <- z$recur==1
early <- !rec & z$death==1 & z$dtime<=z$rtime
late <- !rec & z$death==1 & z$dtime>z$rtime
pt <- ifelse(rec,z$rtime,ifelse(early,z$dtime,z$rtime)); pe <- as.integer(rec|early)
dt <- z$dtime; de <- z$death
separate <- de==0 | ifelse(pe==1,pt<=dt,pt<dt)
common <- pt<=dt & (pe==1 | (de==0 & pt==dt))
stopifnot(sum(rec)==1518, sum(early)==152, sum(late)==43, sum(pe==0)==1312,
          all(pt[late]==z$rtime[late]), all(pe[late]==0), all(separate), sum(!common)==43,
          sum(pt>=tau)==1571, sum(dt>=tau)==2084)
summary <- data.frame(check=c("N","unique_pid","missing_cells","recurrence_branch","composite_death_branch","P_censored","late_deaths_retained","actual_separate_feasible","actual_common_forbidden","P_risk_tau","D_risk_tau","tau_days"),
                    value=c(n,length(unique(z$pid)),sum(is.na(z)),sum(rec),sum(early),sum(pe==0),sum(late),sum(separate),sum(!common),sum(pt>=tau),sum(dt>=tau),tau))
write.csv(summary,args[2],row.names=FALSE)
unlink(folder,recursive=TRUE)
cat("Recorded endpoint construction verified; no survival analysis performed.\n")
