import::from("magrittr", "%>%", "%<>%")
datadir <- file.path(here::here(), "data")
raw_situps <- readr::read_csv(file.path(datadir, "situps 2020-06-2317.17.36.csv"))[1:7]
raw_crunches <- readr::read_csv(file.path(datadir, "crunches 2020-06-2317.13.34.csv"))[1:7]

# I don't really trust the process that generated these data. They have come out dirtier than I would
# really have liked, with loss of precision on the time stamps and irregular sampling. Since these 
# are for demonstration purposes I won't stress too much about it -- I'll simply average over 
# duplicate time stamps and linearly interpolate over gaps, because I'm going to downsample 
# from 1000 Hz to 25 Hz anyway.

regularize_time_sampling <- function(dat) {
  tt <- dat$time
  new_tt <- seq(min(tt), max(tt), by = 0.001)
  apply(dat, 2, function(y) approx(tt, y, new_tt, ties = mean)$y) %>% as.data.frame
  dat$time %<>% round(3)
  dat
}

# Filter out frequencies above 10Hz, more or less
filter_by_brute_force <- function(y, fs = 1000) {
  fmax <- floor(1 + length(y) / 2)
  ff <- seq(0, fs / 2, length.out = fmax)
  ffy <- fft(y)[1:fmax] * pnorm(ff, 11, 1/3, lower.tail = FALSE)
  ffyr <- Conj(rev(ffy[-1]))
  if(length(y) %% 2 == 0) ffyr <- ffyr[-1]
  Re(fft(c(ffy, ffyr), inverse = TRUE))/length(y)
}

# Downsample: filter and decimate (well in this case, quadragintamate)
downsample_from_1000Hz_to_25Hz <- function(dat) {
  dat[2:7] <- lapply(dat[2:7], filter_by_brute_force)
  dat <- as.data.frame(dat)
  dat[seq(1, nrow(dat), by = 40),]
}

(
  .
  %>% regularize_time_sampling
  %>% downsample_from_1000Hz_to_25Hz
) -> preprocess_exercise_data

situps <- preprocess_exercise_data(raw_situps)
crunches <- preprocess_exercise_data(raw_crunches)
saveRDS(situps, "situps.Rds")
saveRDS(crunches, "crunches.Rds")

