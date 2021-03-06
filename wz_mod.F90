!> cleaned
module wz_m
!     Gauss-Labotto and Gauss points
  use kinds, only : DP
  use size_m
  implicit none

  real(DP) :: ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3) &
    ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)

!    Weights

  real(DP) :: WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1) &
    ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2) &
    ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3) &
    ,WAM1(LY1), WAM2(LY2), WAM3(LY3) &
    ,W2AM1(LX1,LY1), W2CM1(LX1,LY1) &
    ,W2AM2(LX2,LY2), W2CM2(LX2,LY2) &
    ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)

end module wz_m
