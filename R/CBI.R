#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' CBI
#'
#' CBI estimates the bias index of a survival dataset attributed to censoring.
#'
#' @param dataset The survival dataset with columns c(time, status, ...)
#' @param time The name of the column corresponding to patient follow-up time.
#' @param status The name of the column corresponding to patient status.
#' @param event A list of labels corresponding to an event in the status column.
#'
#' @return The result is the estimated censored bias index in the interval [0, 1].
#' A threshold must be stablished.
#' @export
#' @references
#' Catalogue of bias collaboration. Bankhead CR, Spencer EA, Nunan D. Information bias.
#'  In: Sackett Catalogue Of Biases 2019.
#'  \href{https://catalogofbias.org/biases/information-bias}{Information bias}
#'
#'  Barrajón E and Barrajón L (2020). Effect of right censoring bias on survival analysis.
#'  \href{https://arxiv.org/abs/2012.08649}{arXiv:2012.08649v1}
#' @examples
#' aml_cbi <- CBI(survival::aml)
#'
#' \dontrun{
#'  ovarian_cbi <- CBI(survival::ovarian) }
#' ovarian_cbi <- CBI(survival::ovarian, time = "futime", status = "fustat")
#'
#' lung_cbi <- CBI(survival::lung, event = 2)
#' veteran_cbi <- CBI(survival::veteran)
#'
CBI <- function(dataset, time = "time", status = "status", event = 1) {
      # try to split in more easier funcitons
      events <- dataset[ dataset[, status] %in% event, time]
      censor <- dataset[!dataset[, status] %in% event, time]
      firstEvent <- min(events)
      lastEvent  <- max(events)
      LT <- sum(censor >= lastEvent)
      censored <- censor[censor >= firstEvent & censor < lastEvent]
      censor2LT <- floor(length(censored) * LT / (LT + length(events)))
      n <- length(censored) - censor2LT
      if(n < 1) return(0)
      atRisk <- head(censored, length(censored) - censor2LT)
      bbi <- sum(atRisk < mean(events)) / length(atRisk)

      adjustSimmetry <- median(censored) / lastEvent
      adjustLT <- length(events) / (LT + censor2LT + length(events))
      bbi ^ (adjustLT * adjustSimmetry)
}



