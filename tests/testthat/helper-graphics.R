skip_if_no_displaylist <- function() {
  has_dl <- tryCatch({
    cap <- grDevices::dev.capabilities()
    isTRUE(cap[["displaylist"]] == "yes")
  }, error = function(e) FALSE)
  if (!has_dl) testthat::skip("Graphics device has no display list; skipping recordPlot tests.")
}
