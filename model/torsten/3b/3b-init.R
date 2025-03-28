# Create Stan initial values
#
# This function must return something that can be passed to the `init` argument
#   of `cmdstanr::sample()`. There are several options; see `?cmdstanr::sample`
#   for details.
#
# `.data` represents the list returned from `make_standata()` for this model.
#   This is provided in case any of your initial values are dependent on some
#   aspect of the data (e.g. the number of rows).
#
# `.args` represents the list of attached arguments that will be passed through to
#   cmdstanr::sample(). This is provided in case any of your initial values are
#   dependent on any of these arguments (e.g. the number of chains).
#
# Note: you _don't_ need to pass anything to either of these arguments, you only
#   use it within the function. `bbr` will pass in the correct objects when it calls
#   `make_init()` under the hood.
#

make_init <- function(.data, .args) {
  function(){
    list(CLHat = exp(rnorm(1, log(45.3),0.25)),
         QHat  = exp(rnorm(1, log(38.3),0.25)),
         V1Hat = exp(rnorm(1, log(7.90),0.25)),
         V2Hat = exp(rnorm(1, log(481),0.25)), 
         KAHat = exp(rnorm(1, log(1),0.25)),
         KA14Hat = exp(rnorm(1, log(5),0.25)),
         KA41Hat = exp(rnorm(1, log(0.1),0.25)),
         nu = max(3.1,rgamma(1,2,0.1)),
         omega = rep(0.5,7) * exp(rnorm(7, 0, 0.25)),
         L = diag(7),
         sigma = exp(rnorm(1,log(0.2),0.25)),
         etaStd = matrix(0L, 7, .data$nSubjects))
  }
}
