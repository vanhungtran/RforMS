#' Common positive- and negative-mode ESI adduct mass shifts
#' @keywords internal
adduct_shifts <- list(
  "[M+H]+"     = 1.007276,
  "[M+Na]+"    = 22.989218,
  "[M+K]+"     = 38.963158,
  "[M+NH4]+"   = 18.033823,
  "[M-H]-"     = -1.007276,
  "[M+Cl]-"    = 34.969402,
  "[M+HCOO]-"  = 44.998201
)

#' Calculate expected m/z for a given neutral mass and adducts
#'
#' @param neutral_mass Numeric exact monoisotopic mass in Daltons
#' @param adducts Character vector of adduct types (e.g., "[M+H]+", "[M+Na]+")
#' @return A named numeric vector of expected m/z values
#' @export
calc_adduct_mz <- function(neutral_mass, adducts = c("[M+H]+", "[M+Na]+")) {
  vapply(adducts, function(add) {
    shift <- adduct_shifts[[add]]
    if (is.null(shift)) stop("Unknown adduct: ", add)
    neutral_mass + shift
  }, numeric(1))
}

#' Convert observed m/z back to putative neutral mass
#'
#' @param mz Observed m/z value
#' @param adduct Assumed adduct type
#' @return Putative monoisotopic neutral mass
#' @export
mz_to_neutral_mass <- function(mz, adduct = "[M+H]+") {
  shift <- adduct_shifts[[adduct]]
  if (is.null(shift)) stop("Unknown adduct: ", adduct)
  mz - shift
}

#' Match observed masses against reference database within ppm tolerance
#'
#' @param observed_masses Numeric vector of observed neutral masses
#' @param reference_db Data frame with 'name' and 'exact_mass' columns
#' @param ppm_tolerance Allowed mass error in parts per million (default: 5 ppm)
#' @return A data frame containing matched hits with calculated ppm errors
#' @export
match_exact_mass <- function(observed_masses, reference_db, ppm_tolerance = 5) {
  matches <- list()

  for (i in seq_along(observed_masses)) {
    obs <- observed_masses[i]
    delta_ppm <- abs(reference_db$exact_mass - obs) / reference_db$exact_mass * 1e6
    hit_idx <- which(delta_ppm <= ppm_tolerance)

    if (length(hit_idx) > 0) {
      matches[[length(matches) + 1]] <- data.frame(
        observed_mass = obs,
        matched_name = reference_db$name[hit_idx],
        reference_mass = reference_db$exact_mass[hit_idx],
        ppm_error = delta_ppm[hit_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(matches) == 0) {
    return(data.frame(
      observed_mass = numeric(), matched_name = character(),
      reference_mass = numeric(), ppm_error = numeric(), stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, matches)
}
