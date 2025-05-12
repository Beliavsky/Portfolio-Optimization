module portfolio_opt_mod
use kind_mod, only: dp
implicit none
private
public :: max_sharpe_long_only, sort_desc, sharpe_ratio, sharpe_ratio_grad

contains

subroutine max_sharpe_long_only(mu, cov, w, sharpe, n, tol, max_iter)
! long-only portfolio that maximizes the ex-ante Sharpe ratio
integer,  intent(in)  :: n               ! number of assets
real(dp), intent(in)  :: mu(n)           ! expected returns
real(dp), intent(in)  :: cov(n,n)        ! covariance matrix
real(dp), intent(out) :: w(n)            ! optimal weights (∑w=1, w≥0)
real(dp), intent(out) :: sharpe          ! optimal Sharpe ratio
real(dp), intent(in),  optional :: tol   ! convergence tolerance
integer,  intent(in),  optional :: max_iter ! maximum iterations

real(dp) :: tol_, ret, var, old_sharpe, alpha, alpha0
integer      :: max_iter_, iter
real(dp), allocatable :: g(:), cw(:), w_new(:)
real(kind=dp), parameter :: alpha_min = 1.0e-20_dp, sharpe_diff = 1.0e-20_dp
tol_      = 1.0e-12_dp           ! default tolerance
if (present(tol)) tol_ = tol
max_iter_ = 10000                   ! default maximum iterations
if (present(max_iter)) max_iter_ = max_iter

allocate(g(n), cw(n), w_new(n))

!---------------- initial equal-weight portfolio --------------------
w = 1.0_dp / n
call project_simplex(w, n)           ! just in case

cw  = matmul(cov, w)
var = dot_product(w, cw)
ret = dot_product(w, mu)
sharpe = ret / sqrt(var)

alpha0 = 0.1_dp                  ! initial step size

!---------------- projected-gradient ascent ------------------------
do iter = 1, max_iter_
   g = mu / sqrt(var) - (ret / var**1.5_dp) * cw  ! Sharpe gradient
   alpha = alpha0
   old_sharpe = sharpe

   do                                  ! backtracking line search
      w_new = w + alpha * g
      call project_simplex(w_new, n)    ! project to simplex (long-only, fully-invested)

      cw  = matmul(cov, w_new)
      var = dot_product(w_new, cw)
      ret = dot_product(w_new, mu)
      sharpe = ret / sqrt(var)

      if (sharpe > old_sharpe + sharpe_diff .or. alpha < alpha_min) exit
      alpha = alpha * 0.5_dp
   end do

   if (abs(sharpe - old_sharpe) < tol_) exit
   w = w_new
end do
print*,"in max_sharpe_long_only, exiting loop, iter, alpha =",iter,alpha ! debug
end subroutine max_sharpe_long_only

subroutine project_simplex(v, n)
! projects vector v onto the probability simplex {x | x≥0, Σx=1}
integer , intent(in)     :: n
real(dp), intent(in out) :: v(n)
real(dp), allocatable    :: u(:)
real(dp) :: theta, cumsum
integer :: k

u = v
call sort_desc(u)

cumsum = 0.0_dp
do k = 1, n
   cumsum = cumsum + u(k)
   theta = (cumsum - 1.0_dp) / k
   if (k == n) exit
   if (u(k+1) - theta <= 0.0_dp) exit
end do

v = max(v - theta, 0.0_dp)
end subroutine project_simplex

subroutine sort_desc(a)
! insertion sort in descending order (small n, so O(n^2) is fine)
real(dp), intent(inout) :: a(:)
integer :: i, j, n
real(dp) :: key
n = size(a)
do i = 2, n
   key = a(i)
   j = i - 1
   do while (j >= 1)
      if (a(j) >= key) exit
      a(j+1) = a(j)
      j = j - 1
   end do
   a(j+1) = key
end do
end subroutine sort_desc

function sharpe_ratio(w, mu, cov) result(sharpe)
! computes the ex-ante Sharpe ratio for weights w given mu and cov
real(dp), intent(in) :: w(:)         ! portfolio weights
real(dp), intent(in) :: mu(:)        ! expected returns
real(dp), intent(in) :: cov(size(w),size(w)) ! covariance matrix
real(dp) :: sharpe, ret, var
real(dp), allocatable :: cw(:)
cw  = matmul(cov, w)
var = dot_product(w, cw)
ret = dot_product(w, mu)
sharpe = ret / sqrt(var)
end function sharpe_ratio

subroutine sharpe_ratio_grad(w, mu, cov, sharpe, grad)
! computes the Sharpe ratio and its gradient wrt w
real(dp), intent(in)  :: w(:)                       ! weights
real(dp), intent(in)  :: mu(:)                      ! expected returns
real(dp), intent(in)  :: cov(size(w),size(w))       ! covariance matrix
real(dp), intent(out) :: sharpe                     ! Sharpe ratio
real(dp), intent(out) :: grad(size(w))              ! gradient dS/dw
real(dp) :: ret, var
real(dp), allocatable :: cw(:)
! compute ret and var
cw = matmul(cov, w)
var = dot_product(w, cw)
ret = dot_product(w, mu)
sharpe = ret / sqrt(var)
grad = mu / sqrt(var) - (ret / var**1.5_dp) * cw
end subroutine sharpe_ratio_grad

end module portfolio_opt_mod
