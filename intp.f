!* Author: Michele Trabucchi

!* THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!* SOFTWARE.

!* Interpolation program
      program intp
      implicit none
      integer tab_size,       !: number of rows in the table
     `        expected_narg,  !: expected number of arguments
     `        uTAB,uIN,uOUT   !: units for table, input, output files
      parameter(tab_size=324,expected_narg=3,uTAB=1,uIN=10,uOUT=11)
      !
      integer narg,     !: number of received command-line arguments
     `        n_input,  !: number of rows in the input file
     `        i,j
      character(1) mode !: mode (use hydrostatic/hydrodynamic values)
      character(256) line,            !: line read from the input file
     `               input_filename,  !: name of the input file
     `               output_filename  !: name of the output file
      logical input_exists,output_exists
      !
      ! Grid structure and table
      !;    X  -> hydrogen mass fraction
      !;    Z  -> metallicity by mass fraction
      !;    M  -> stellar mass in solar units
      !;    AM -> mixing length parameter (α_ML)
      !;    AN -> turbulent viscosity parameter (α_ν)
      ! Number of nodes in each dimension
      integer nX,nZ,nAM,nM,nAN
      parameter(nX=3,nZ=3,nAM=3,nM=2,nAN=6)
      ! Values at the grid nodes in each dimension
      real*8 kX(nX),kZ(nZ),kAM(nAM),kM(nM),kAN(nAN) !: grid nodes
      parameter(kX=(/0.6d0,0.7d0,0.8d0/))
      parameter(kZ=(/0.002d0,0.006d0,0.010d0/))
      parameter(kAM=(/1.0d0,1.5d0,2.0d0/))
      parameter(kM=(/1.0d0,2.6d0/))
      parameter(kAN=(/0.0d0,0.01d0,0.02d0,0.05d0,0.10d0,0.20d0/))
      ! Bracketing parameters for evenly-spaced nodes
      real*8 gdX,gdZ,gdM !: grid spacing by value
      real*8 gnX,gnZ,gnM !: grid spacing by number of nodes
      ! Table of values for interpolation
      real*8 gX(tab_size),gZ(tab_size),gM(tab_size),
     `       gAM(tab_size),gAN(tab_size),
     `       gLs(tab_size),gTs(tab_size),
     `       gLd(tab_size),gTd(tab_size),
     `       gndom(tab_size),gP(tab_size)
      !
      ! Variables used in the interpolation
      real*8 p1,p2,p3     !: Lagrange basis polynomials for quadratic interpolation
      real*8 w1AN,w1M,w1Z,w1X,w2AN,w2M,w2Z,w2X !: linear interpolation weights
      integer*1 fX,fZ,fM,fAN,fT !: flags for linear interpolation
      !
      real*8 tT(3)        !: log(Te) values corresponding to the three α_ML nodes at fixed X,Z,M,α_ν, used to compute Lagrange basis polynomials
      real*8 tL(3),tP(3)  !: log(L) and log(P) values corresponding to the three α_ML nodes at fixed X,Z,M,α_ν, in quadratic interpolation formula
      real*8 aL(2),aP(2)  !: stores the results of quadratic interpolation in log(Te) at fixed X,Z,M for two values of α_ν, used to interpolate linearly in the latter
      real*8 mL(2),mP(2)  !: stores the results of linear interpolation in α_ν at fixed X,Z for two values of M, used to interpolate linearly in the latter
      real*8 zL(2),zP(2)  !: stores the results of linear interpolation in M at fixed X for two values of Z, used to interpolate linearly in the latter
      real*8 xL(2),xP(2)  !: stores the results of linear interpolation in Z for two values of X, used to interpolate linearly in the latter
      real*8 logL,logP    !: stores final value of log(L) and log(P) interpolated in X,Z,M,α_ν,log(Te)
      !
      real*8 X,Z,M,AN,lTe !: input values of X,Z,M,α_ν,log(Te)
      integer jX,jZ,jM,jAM,jAN !: running indices for the grid nodes
      integer jX0,jZ0,jM0,jAN0 !: indices of the lower bracketing node
      !
  500 format(f3.1,f6.3,f4.1,f4.1,f5.2,f6.3,f5.0,f6.3,f5.0,f3.0,f7.1)
  600 format("Missing arguments!")
  601 format("Bad arguments!")
  602 format("Usage: <'s' or 'd'> <input_filename>")
  603 format('File not found: ',a)
  604 format('Output file already exists: ',a)
  605 format('Reading input from file: ',a)
  606 format('Output written to file: ',a)
  699 format(a,f7.3,'%')
  700 format(a,1x,2es24.15e3,5i2)
      !
      ! Read table of values from file
      ! (by default use unit 1 file "fort.1" in the current directory)
      read(uTAB,fmt='(a)')
      do i=1,tab_size
        read(uTAB,fmt=500) gX(i),gZ(i),gAM(i),gM(i),gAN(i),
     `                     gLs(i),gTs(i),gLd(i),gTd(i),gndom(i),gP(i)
        gTs(i)=log10(gTs(i))
        gTd(i)=log10(gTd(i))
        gP(i)=log10(gP(i))
      end do
      ! Compute distances for evenly-spaced nodes
      gnX=real(nX-1,8)
      gnZ=real(nZ-1,8)
      gnM=real(nM-1,8)
      gdX=kX(nX)-kX(1)
      gdZ=kZ(nZ)-kZ(1)
      gdM=kM(nM)-kM(1)
      !
! ======================================================================
      ! Read command-line arguments
      narg = iargc()
      if (narg.ne.expected_narg) then
        write(6,600)
        write(6,602)
        stop
      end if
      call getarg(1,mode)
      if (mode.ne.'s'.and.mode.ne.'d') then
        write(6,601)
        write(6,602)
        stop
      end if
      call getarg(2,input_filename)
      inquire(file=trim(input_filename),exist=input_exists)
      if (.not.input_exists) then
        write(6,603) trim(input_filename)
        stop
      end if
      call getarg(3,output_filename)
      inquire(file=trim(output_filename),exist=output_exists)
      if (output_exists) then
        write(6,604) trim(output_filename)
        stop
      end if
      write(6,605) trim(input_filename)
      open(unit=uIN,file=trim(input_filename),
     `     form='formatted',status='old')
      n_input=0
      do
        read(unit=uIN,fmt=*,end=900)
        n_input=n_input+1
      enddo
  900 close(unit=uIN)
! ======================================================================
      ! Open input/output files
      open(unit=uIN,file=trim(input_filename),
     `     form='formatted',status='old')
      open(unit=uOUT,file=trim(output_filename),
     `     form='formatted',status='new')
      do i=1,n_input
        ! Read input values at which to interpolate
        read(unit=uIN,fmt='(a)') line
        read(line,fmt=*) X, Z, M, AN, lTe
        ! Attempt to infer if Teff is in linear scale, and convert it
        if (lTe.ge.100.d0) then
          lTe=log10(lTe)
        endif
        ! Bracket nodes in X, Z, M (regular spacing)
        jX0=max(min(int(floor(gnX*(X-kX(1))/gdX))+1,2),1)
        jZ0=max(min(int(floor(gnZ*(Z-kZ(1))/gdZ))+1,2),1)
        jM0=1
        !; jM0=max(min(int(floor(gnM*(M-kM(1))/gdM))+1,2),1) ! formally correct but redundant
        ! Bracket nodes in α_ν (irregular spacing)
        do jAN=2,nAN
          if (AN.le.kAN(jAN)) then
            jAN0=jAN-1
            exit
          endif
        end do
        ! Interpolation flags
        if (X.lt.kX(1)) then
          fX=1
        else if (X.gt.kX(nX)) then
          fX=2
        else
          fX=0
        endif
        if (Z.lt.kZ(1)) then
          fZ=1
        else if (Z.gt.kZ(nZ)) then
          fZ=2
        else
          fZ=0
        endif
        if (M.lt.kM(1)) then
          fM=1
        else if (M.gt.kM(nM)) then
          fM=2
        else
          fM=0
        endif
        if (AN.lt.kAN(1)) then
          fAN=1
        else if (AN.gt.kAN(nAN)) then
          fAN=2
        else
          fAN=0
        endif
        ! Linear interpolation weights
        w1X=(kX(jX0+1)-X)/(kX(jX0+1)-kX(jX0))
        w2X=1.d0-w1X
        w1Z=(kZ(jZ0+1)-Z)/(kZ(jZ0+1)-kZ(jZ0))
        w2Z=1.d0-w1Z
        w1M=(kM(jM0+1)-M)/(kM(jM0+1)-kM(jM0))
        w2M=1.d0-w1M
        w1AN=(kAN(jAN0+1)-AN)/(kAN(jAN0+1)-kAN(jAN0))
        w2AN=1.d0-w1AN
        ! Multi-dimensional interpolation
        do jX=1,2
          do jZ=1,2
            do jM=1,2
              do jAN=1,2
                do jAM=1,3
                  ! Exploit table structure to find the three rows that
                  ! list the values of log(L), log(Te) for the required
                  ! combination of X,Z,α_ML,M,α_ν
                  !; Complete form: j=1+108*(jX0+jX-2)+36*(jZ0+jZ-2)+12*(jAM-1)+6*(jM0+jM-2)+(jAN0+jAN-2)
                  j=108*(jX0+jX)+36*(jZ0+jZ)+12*(jAM)
     `               +6*(jM0+jM)+jAN0+jAN-313
                  ! Gather the three values of log(L), log(P), log(Te)
                  ! that correspond to the three values of α_ML
                  tP(jAM)=gP(j)
                  if (mode.eq.'s') then
                    ! Select the hydrostatic values
                    tL(jAM)=gLs(j)
                    tT(jAM)=gTs(j)
                  else
                    ! Select the hydrodynamic time-average values
                    tL(jAM)=gLd(j)
                    tT(jAM)=gTd(j)
                  endif
                enddo
                ! Flag for quadratic interpolation
                if (lTe.lt.tT(1)) then
                  fT=1
                else if (lTe.gt.tT(3)) then
                  fT=2
                else
                  fT=0
                endif
                ! Compute Lagrange basis polynomials for quadratic interpolation
                p1=(lTe-tT(2))*(lTe-tT(3))/(tT(1)-tT(2))/(tT(1)-tT(3))
                p2=(lTe-tT(1))*(lTe-tT(3))/(tT(2)-tT(1))/(tT(2)-tT(3))
                p3=(lTe-tT(1))*(lTe-tT(2))/(tT(3)-tT(1))/(tT(3)-tT(2))
                ! Interpolate quadratically log(L) and log(P) in log(Te)
                aL(jAN)=p1*tL(1)+p2*tL(2)+p3*tL(3)
                aP(jAN)=p1*tP(1)+p2*tP(2)+p3*tP(3)
                ! Two values are obtained for each node of α_ν that
                ! bracket the input α_ν, they are stored in aL for
                ! linear interpolation in α_ν
              enddo
              ! Perform linear interpolation of log(L) in α_ν
              mL(jM)=w1AN*aL(1)+w2AN*aL(2)
              mP(jM)=w1AN*aP(1)+w2AN*aP(2)
              ! Two values are obtained for each node of M that
              ! bracket the input M, they are stored in mL for
              ! linear interpolation in M
            enddo
            ! Perform linear interpolation of log(L) in M
            zL(jZ)=w1M*mL(1)+w2M*mL(2)
            zP(jZ)=w1M*mP(1)+w2M*mP(2)
            ! Two values are obtained for each node of Z that
            ! bracket the input Z, they are stored in zL for
            ! linear interpolation in Z
          enddo
          ! Perform linear interpolation of log(L) in Z
          xL(jX)=w1Z*zL(1)+w2Z*zL(2)
          xP(jX)=w1Z*zP(1)+w2Z*zP(2)
          ! Two values are obtained for each node of X that
          ! bracket the input X, they are stored in xL for
          ! linear interpolation in X
        enddo
        ! Perform linear interpolation of log(L) in X
        logL=w1X*xL(1)+w2X*xL(2)
        logP=w1X*xP(1)+w2X*xP(2)
        ! Write interpolated value to output file along with input
        write(uOUT,700) trim(line),logL,logP,fX,fZ,fM,fAN,fT
        call flush(6)
        write(6,699,advance='no') achar(13),real(i)/real(n_input)*100.0
      enddo
      write(6,*) ' DONE!'
      close(unit=uIN)
      close(unit=uOUT)
      write(6,606) trim(output_filename)
      !
      stop
      end program intp