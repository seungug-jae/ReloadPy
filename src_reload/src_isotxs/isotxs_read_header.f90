   subroutine isotxs_read_header(fileunit,user_isotxs,iopt)

!  Feb. 2011 : Initial production version, C. H. Lee
!              Modified from the UNIC routines

!  iopt  = 0 : read the first 2 records and do backspaces
!          1 : read all 4 records
!          2 : read all 4 records and do backspaces

   use m_kind
   use m_isotxs
   
   implicit none
   
   integer, intent(IN) :: fileunit,iopt
   type (type_isotxs)  :: user_isotxs
!  local
   integer i,ig,ngroup,niso,ichist

!  card type 0
   read(fileunit) user_isotxs%hname,(user_isotxs%huse(i),i=1,2),user_isotxs%ivers

!  card type 1: file control (1d record)
   read(fileunit)  user_isotxs%ngroup,user_isotxs%niso,user_isotxs%maxup,user_isotxs%maxdn,user_isotxs%maxord, &
                   user_isotxs%ichist,user_isotxs%nscmax,user_isotxs%nsblok
   if(iopt == 0) then
     do i=1,2
       backspace fileunit
     enddo
     return
   endif

   ngroup=user_isotxs%ngroup
   niso=user_isotxs%niso
   ichist=user_isotxs%ichist

!  card type 2: file data (2d record)
   if(user_isotxs%ichist == 1) then
     allocate(user_isotxs%hisonm(niso),user_isotxs%filechi(ngroup,1), &
              user_isotxs%vel(ngroup),user_isotxs%energy(ngroup+1),user_isotxs%loca(niso))
     read(fileunit)  (user_isotxs%hsetid(i),i=1,12),(user_isotxs%hisonm(i),i=1,user_isotxs%niso),                     &
                     (user_isotxs%filechi(i,1),i=1,user_isotxs%ngroup),(user_isotxs%vel(i),i=1,user_isotxs%ngroup), &
                     (user_isotxs%energy(i),i=1,user_isotxs%ngroup+1),(user_isotxs%loca(i),i=1,user_isotxs%niso)
   else if(user_isotxs%ichist == 0) then
     allocate(user_isotxs%hisonm(niso),user_isotxs%vel(ngroup), &
              user_isotxs%energy(ngroup+1),user_isotxs%loca(niso))
     read(fileunit)  (user_isotxs%hsetid(i),i=1,12),(user_isotxs%hisonm(i),i=1,user_isotxs%niso),                     &
                     (user_isotxs%vel(i),i=1,user_isotxs%ngroup),                                                   &
                     (user_isotxs%energy(i),i=1,user_isotxs%ngroup+1),(user_isotxs%loca(i),i=1,user_isotxs%niso)
   endif

!  card type 3: file-wide chi data (3d record)
   if(user_isotxs%ichist > 1) then
     allocate(user_isotxs%filechi(ngroup,ichist))
     read(fileunit) ((user_isotxs%filechi(ig,i),ig=1,user_isotxs%ngroup),i=1,user_isotxs%ichist)
   endif

   if(iopt == 2) then
     do i=1,4
       backspace fileunit
     enddo
     return
   endif

!  deallocate

   end subroutine
