program xportfolio_opt
  use kind_mod         , only: dp
  use portfolio_opt_mod, only: max_sharpe_long_only, sort_desc, &
                               sharpe_ratio
  implicit none

  integer, parameter :: n = 100        ! <-- change n here
  real(dp)           :: mu(n)          ! simulated expected returns
  real(dp)           :: cov(n,n)       ! simulated covariance matrix
  real(dp)           :: w(n)           ! optimal weights
  real(dp)           :: sharpe         ! optimal Sharpe ratio
  real(dp)           :: A(n,n)
  real(dp), parameter :: low = 0.08_dp, high = 0.10_dp

  ! initialize random-number generator
  call random_seed()

  call random_number(mu)
  mu = low + (high - low) * mu
  call sort_desc(mu)

  ! simulate covariance = A^T * A
  call random_number(A)
  cov = matmul(transpose(A), A) / real(n, dp)

  ! perform long-only maximum-Sharpe optimization
  call max_sharpe_long_only(mu, cov, w, sharpe, n)

  print "(/,a,/,*(f8.4))", "Simulated expected returns (mu):", mu
  print "(/,a,/,*(f8.4))", "Optimal long-only weights (sum = 1):", w
  print "(/,a, f8.4)", "Maximum Sharpe ratio: ", sharpe
  print "(/,a, f8.4)", "Check Sharpe ratio: ", sharpe_ratio(w, mu, cov)
  print*,"#nonzero weights:", count(w > 0.0_dp)

end program xportfolio_opt
