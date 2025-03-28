make_standata <- function(.dir) {
  ## get data file
  xdata <- readr::read_csv(here::here("data", "derived", "nonmem_data_v02.csv "))

  xdata <- xdata %>%
    filter(OUT!=1, CMT!=22, CMT!=33)
  
  #' Format data for Stan
  nt   <- nrow(xdata)
  iObs <- with(xdata, (1:nrow(xdata))[EVID == 0])
  nObs <- length(iObs)
  xsub <- subset(xdata, !duplicated(ID))
  nSubjects <- length(xsub$ID)
  start <- (1:nt)[!duplicated(xdata$ID)]
  end <- c(start[-1] - 1, nt)

  data <- with(xdata,
                    list(nt = nt,
                         nObs = nObs,
                         nSubjects = nSubjects,
                         iObs = iObs,
                         start = start,
                         end = end,
                         time = TIMEh,
                         cObs = DV[iObs],
                         amt =  AMT,
                         rate = RATE,
                         cmt = CMT,
                         evid = EVID,
                         ii = 0*CMT,
                         addl = 0*CMT,
                         ss = 0*CMT,
                         runestimation=1))
}