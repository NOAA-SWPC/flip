!date 20150810
!written by Tim Fuller Rowell
!purpose: initialization of tiros routine
      subroutine tiros_init()
      use tirosdata
      CHARACTER (LEN=*), parameter :: filename3='ionprof'
      CHARACTER (LEN=*), parameter :: filename7='tiros_spectra'
      CHARACTER (LEN=100) :: string_dum
      INTEGER :: istat, i
      INTEGER, parameter :: UNIT3=103
      INTEGER, parameter :: UNIT7=107
      open(UNIT=UNIT3,file=trim(filename3),status='old',
     &form='formatted', iostat=istat)
      IF ( istat /= 0 ) THEN
         WRITE( UNIT=6, FMT=*)'ERROR OPENING FILE',filename3
         STOP
      END IF
      open(UNIT=UNIT7,file=trim(filename7),status='old',
     &form='formatted', iostat=istat)
      IF ( istat /= 0 ) THEN
         WRITE( UNIT=6, FMT=*)'ERROR OPENING FILE',filename7
         STOP
      END IF

99001 FORMAT (1x,6E13.6)
         READ (103,99001) emaps
         READ (103,99001) cmaps
         read(UNIT=UNIT7,fmt=*) string_dum
         read(UNIT=UNIT7,fmt=*) string_dum
         read(107,fmt=*)
         do iband=1,21
         read(107,fmt=*) string_dum
         read(107,fmt=*)
         read(UNIT=UNIT7,fmt='(1X,5e10.4)')(djspectra(iflux,iband),
     &iflux=1,15)
         read(107,fmt=*)
         enddo
      CLOSE(UNIT=UNIT3)
      CLOSE(UNIT=UNIT7)
         return
      end subroutine tiros_init
       
