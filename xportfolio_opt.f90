program xportfolio_opt
use kind_mod, only: dp
use portfolio_opt_mod , only: max_sharpe_long_only
implicit none

integer , parameter :: n = 4
real(dp) :: mu(n)            ! expected returns
real(dp) :: cov(n,n)         ! covariance matrix
real(dp) :: w(n)             ! optimal weights
real(dp) :: sharpe           ! optimal Sharpe ratio

! example inputs -----------------------------------------------------
mu  = [0.12_dp, 0.10_dp, 0.07_dp, 0.03_dp]

cov = reshape([0.10_dp, 0.01_dp, 0.02_dp, 0.015_dp, &
               0.01_dp, 0.08_dp, 0.018_dp, 0.012_dp, &
               0.02_dp, 0.018_dp, 0.07_dp, 0.010_dp, &
               0.015_dp, 0.012_dp, 0.010_dp, 0.05_dp ], &
              [n, n])

! long-only maximum-Sharpe optimization ------------------------------
call max_sharpe_long_only(mu, cov, w, sharpe, n)

! results ------------------------------------------------------------
print *, "optimal long-only weights (sum = 1):"
print "(*(f8.4))", w
print "(a,f8.4)", "maximum Sharpe ratio = ", sharpe

end program xportfolio_opt
