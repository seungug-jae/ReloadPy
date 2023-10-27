   module m_kind

!  Feb. 2011 : Initial production version, C. H. Lee

   integer, parameter :: sp = selected_real_kind(p =  6, r =  37)               ! single-precision 
   integer, parameter :: dp = selected_real_kind(p = 13, r = 200)               ! double-precision 
!  integer, parameter :: qp = selected_real_kind(p = 26, r = 300)               ! quadruple-precision 

   integer, parameter :: singel    = sp                            
   integer, parameter :: doubel    = dp                            
!  integer, parameter :: quadrupel = qp                            

   integer, parameter :: highp     = dp                            
   integer, parameter :: doubles   = dp
   integer, parameter :: reals     = sp                           

   end module
