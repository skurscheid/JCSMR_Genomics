# In RepeatReads.R:
#' @docType package
#' 
#' RepeatRead class description
#' @import methods
#' @export RepeatRead
#' @exportClass RepeatRead

RepeatRead <- setRefClass("RepeatRead",
  fields = list(repClass = "numeric",
                repFamily = "numeric"),
  methods = list(
    unique = function(x) {
      if (sum(repClass) == 1) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    },
    rep.Class = function(x) {
      names(repClass[which(repClass == 1)])
    },
    rep.Family = function(x) {
      names(repFamily[which(repFamily == 1)])
    }
  ))
