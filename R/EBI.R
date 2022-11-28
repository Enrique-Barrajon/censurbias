#' EBI
#'
#' EBI estimates the bias index of a survival dataset attributed to censoring.
#'
#' @param dataset The survival dataset with columns c(time, status, ...)
#' @param time The name of the column corresponding to patient follow-up time.
#' @param status The name of the column corresponding to patient status.
#' @param event A list of labels corresponding to an event in the status column.
#'
#' @return The result is the estimated censored bias index in the interval [0, 1].
#' A value over 0.5 suggests right-censoring information bias.
#' @export
#' @references
#' Catalogue of bias collaboration. Bankhead CR, Spencer EA, Nunan D. Information bias.
#'  In: Sackett Catalogue Of Biases 2019.
#'  \href{https://catalogofbias.org/biases/information-bias}{Information bias}
#'
#'  Barrajón E and Barrajón L (2020). Effect of right censoring bias on survival analysis.
#'  \href{https://arxiv.org/abs/2012.08649}{arXiv:2012.08649v1}
#' @examples
#' aml_cbi <- EBI(survival::aml)
#'
#' \dontrun{
#'  ovarian_cbi <- EBI(survival::ovarian) }
#' ovarian_cbi <- EBI(survival::ovarian, time = "futime", status = "fustat")
#'
#' lung_cbi <- EBI(survival::lung, event = 2)
#' veteran_cbi <- EBI(survival::veteran)
#'
EBI <- function(dataset, time = "time", status = "status", event = 1, threshold = .96) {
      centering <- log(0.5) / log(threshold)
      CBI(dataset, time, status, event) ^ centering
}


