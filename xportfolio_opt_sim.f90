program xportfolio_opt
  use kind_mod         , only: dp
  use portfolio_opt_mod, only: max_sharpe_long_only, sort_desc, &
                               sharpe_ratio, sharpe_ratio_grad
  implicit none

  integer, parameter :: n = 1000       ! <-- change n here
  real(dp)           :: mu(n)          ! simulated expected returns
  real(dp)           :: cov(n,n)       ! simulated covariance matrix
  real(dp)           :: w(n)           ! optimal weights
  real(dp)           :: sharpe         ! optimal Sharpe ratio
  real(dp)           :: sharpe_max_ran
  real(dp)           :: A(n,n)
  real(dp)           :: grad(n)        ! gradient of Sharpe ratio
  real(dp), parameter :: low = 0.01_dp, high = 0.10_dp
  integer, parameter :: ntry_ran_w = 10**4
  integer :: i
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
  print "(a, f8.4)", "Check Sharpe ratio: ", sharpe_ratio(w, mu, cov)
  call sharpe_ratio_grad(w, mu, cov, sharpe, grad)
  print "(/,a,/,*(f8.4))", "Gradient of Sharpe ratio:", grad
  print "(/,a,i0)", "#nonzero weights: ", count(w > 0.0_dp)
  w = 1.0_dp/n
  call sharpe_ratio_grad(w, mu, cov, sharpe, grad)
  print "(/,a, f8.4)", "Sharpe ratio for equal weights: ", sharpe
  print "(/,a,/,*(f8.4))", "Gradient of Sharpe ratio:", grad
  if (ntry_ran_w > 0) then
     sharpe_max_ran = -huge(sharpe_max_ran)
     do i=1,ntry_ran_w
        call random_number(w)
        sharpe = sharpe_ratio(w, mu, cov)
        sharpe_max_ran = max(sharpe, sharpe_max_ran)
     end do
     print "(/,'# random sets of weights: ',i0)", ntry_ran_w
     print "(a, f8.4)", "maximum random Sharpe:", sharpe_max_ran
  end if
end program xportfolio_opt
