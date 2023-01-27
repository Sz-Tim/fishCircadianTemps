

phi_to_ZT <- function(phi) {
  phi <- phi + 2*pi*(phi < -pi) - 2*pi*(phi > pi)
  phi <- (-phi+pi)*12/pi-12
  phi <- phi + 24*(phi < 0)
}


calc_prefTemp <- function(M, A, phi, ZT) {
  M + exp(A) * cos(3.141593*ZT/12 + phi)
}
