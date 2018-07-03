#' Generate a stimulus boxcar function.
#' 
#' Generates a stimulus boxcar vector for the specified time duration
#' and microtime resolution based on the user-defined onsets and durations.
#' 
#' @export
#' @importFrom dplyr if_else funs mutate n rename_at tibble %>%
#' @importFrom purrr map map2 reduce rerun
#' 
#' @param totaltime Total time of the design in seconds.
#' @param onsets Vector representing the onsets of the stimulus in seconds.
#' @param durations Vector representing the durations of the stimulus in seconds.
#' @param cond.names Optional names for the conditions.
#' @param accuracy Microtime resolution in seconds.
#' @return A tibble in microtime resolution specifying the stimulus boxcar
#' functions, with one column per condition, in 1-0 coding.
#' 
#' @details If duration is a single number, it is assumed that all stimulus
#' onsets have the same duration.
#' 
#' @seealso \code{\link{specifydesign}}
#' 
#' @examples
#' total <- 100
#' os <- c(1, 21, 41, 61, 81)
#' d <- 10
#' out <- stimfunction(total, os, d, accuracy = 0.1)
#' de\dontshow{rm(total,os,d,out)}
#' 
#' @keywords low-level

stimfunction <- function(totaltime, onsets, durations, cond.names = NULL, accuracy){
  # TODO: allow stimfunction to output an object that is the input to specifydesign
  # with multiple onset conditions
  
  # these two conditionals should stop when c() or a list with all entries NULL are input for onsets
  if (length(onsets) < 1) stop("Must specify more than zero onsets.")
  if (sum(lengths(onsets)) < 1) stop("Must specify more than zero onsets.")
  
  if (!is.list(onsets)) onsets <- list(onsets) # if input as num vector, get into 1 list field per condition format
  
  if (length(unlist(durations)) == 1) {
    # in this case, all onsets/all conditions same duration
    durations <- map(onsets, ~rep(unlist(durations), length(.)))
  } else if (!is.list(durations)) {
    # each duration corresponds to an onset condition
    if (length(durations) == length(onsets)) {
      durations <- map2(durations, onsets, ~rep(.x, length(.y)))
    } else if (length(onsets) == 1 & length(durations) == lengths(onsets)) {
      # need this special case now that onsets are always coerced to list
      # if one onset condition, and a duration is specified individually for each onset
      durations <- list(durations)
    }
  } else if (is.list(durations) & length(durations) == length(onsets)) {
    # each list-field of duration corresponds to an onset condition
    if (all(lengths(durations) == 1)) {
      # if every onset within condition has the same onset:
      durations <- map2(durations, onsets, ~rep(.x, length(.y)))
    } else if (all(lengths(durations) == lengths(onsets))) {
      # do nothing
    } else {
      stop("Mismatch between onsets and durations.")
    }
  } else {
    stop("Mismatch between onsets and durations.")
  }
  
  if (max(unlist(onsets)) > totaltime) stop("Mismatch between onsets and totaltime.")
  
  if (is.null(cond.names)) cond.names <- paste0("C", 1:length(onsets))
  
  if (length(cond.names) != length(onsets)) stop("Mismatch between # onset conditions & # condition names.")

  
  # There's an under/overflow issue such that we can't rely on using microtime itself as the index column
  out <- tibble(cond.names = cond.names,
                # will calculate onsets/durations scaled by accuracy to get everything into integers
                params = map2(onsets, durations, ~tibble(onset = as.integer(.x/accuracy), duration = as.integer(.y/accuracy)) %>%
                                mutate(micro = map2(onset, duration, ~seq(.x, .x + .y, 1L))))) %>%
    # my understanding is that the index starts after 0
    mutate(full = rerun(n(), tibble(microtime = seq(1L, totaltime/accuracy, 1L))),
           full = map2(full, params, ~.x %>%
                         mutate(s = if_else(microtime %in% unlist(.y$micro), 1L, 0L),
                                # need to get it back out of integers into real seconds
                                microtime = microtime * accuracy)),
           # must de-duplicate names BEFORE the full_join call for it to work with >2 conditions
           full = map2(full, cond.names, ~rename_at(.x, "s", funs(paste(., .y, sep = ".")))))
  
  return(reduce(out$full, full_join, by = "microtime"))
}

