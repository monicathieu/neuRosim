context("Stimulus function generation")
library(neuRosim2)

test_that("stimfunction outputs tibble", {
  expect_s3_class(stimfunction(totaltime = 100,
                               onsets = 1,
                               durations = 10,
                               accuracy = 0.1),
                  "tbl_df")
  expect_s3_class(stimfunction(totaltime = 100,
                               onsets = seq(1, 91, 10),
                               durations = 5,
                               accuracy = 0.1),
                  "tbl_df")
  expect_s3_class(stimfunction(totaltime = 100,
                               onsets = list(seq(1, 81, 20), seq(11, 91, 20)),
                               durations = c(5, 4), accuracy = 0.1),
                  "tbl_df")
})

test_that(">= 1 onset time specified works", {
  stim <- stimfunction(totaltime = 100, onsets = 1, durations = 10, accuracy = 0.1)
  expect_equal((10/0.1)+1, sum(stim$s.C1))
  stim2 <- stimfunction(totaltime = 100, onsets = 0.1, durations = 10, accuracy = 0.1)
  expect_equal((10/0.1)+1, sum(stim$s.C1))
})

test_that("one condition specified un-listed works", {
  stim <- stimfunction(totaltime = 100, onsets = seq(1, 81, 20), durations = 10, accuracy = 0.1)
  expect_equal((10*5/0.1)+(1*5), sum(stim$s.C1))
  stim2 <- stimfunction(totaltime = 100, onsets = seq(1, 81, 20), durations = c(10, 5, 5, 5, 10), accuracy = 0.1)
  expect_equal(((10*2+5*3)/0.1)+(1*5), sum(stim2$s.C1))
})

test_that("multiple conditions specified listed works", {
  stim <- stimfunction(totaltime = 100, onsets = list(1, 21, 41), durations = 10, accuracy = 0.1)
  expect_equal((10*1/0.1)+(1*1), sum(stim$s.C1))
  expect_equal((10*1/0.1)+(1*1), sum(stim$s.C2))
  expect_equal((10*1/0.1)+(1*1), sum(stim$s.C3))
  stim2 <- stimfunction(totaltime = 100, onsets = list(seq(1, 81, 20), seq(11, 91, 20)), durations = list(5, 4), accuracy = 0.1)
  expect_equal(((5*5)/0.1)+(1*5), sum(stim2$s.C1))
  expect_equal(((4*5)/0.1)+(1*5), sum(stim2$s.C2))
})

test_that("0 onsets specified fails", {
  expect_error(stimfunction(totaltime = 100, onsets = c(), durations = 10, accuracy = 0.1))
  expect_error(stimfunction(totaltime = 100, onsets = vector("list", 3), durations = 10, accuracy = 0.1))
})

test_that("duration 0 returns stick-ish function", {
  stim <- stimfunction(totaltime = 100, onsets = seq(1, 91, 10), durations = 0, accuracy = 0.1)
  expect_equal(10, sum(stim$s.C1))
})

expect_in <- function (object, expected_in) {
  expect_true(object %in% expected_in,
              info = paste(object, "was not found in", tail(substitute(expected_in), 1)))
}

test_that("condition name spec behaves as planned", {
  stim <- stimfunction(totaltime = 100, onsets = list(1, 21, 41), durations = 10, cond.names = c("first", "second", "third"), accuracy = 0.1)
  expect_in("s.first", names(stim))
  expect_in("s.second", names(stim))
  expect_in("s.third", names(stim))
  expect_error(stimfunction(totaltime = 100, onsets = list(1, 21, 41), durations = 10, cond.names = c("first", "second"), accuracy = 0.1))
})
