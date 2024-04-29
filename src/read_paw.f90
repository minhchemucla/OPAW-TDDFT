subroutine readpaw(it)
    use param_mod, only : rpad, rpad_r,toll=> tollsij
    use paw_mod
    use main_mod
    use atom_mod, only : p => pawinfo, at => atominfo
    implicit none

    integer :: it
    !for reading
    character(len=150) :: line
    character(len= 25) :: s1="atom"
    character(len= 25) :: s2="<valence_states>"
    character(len= 25) :: s3="</valence_states>"
    character(len= 25) :: s4="radial_grid"
    character(len= 25) :: s5="<projector_function"
    character(len= 25) :: s6="<ae_partial_wave"
    character(len= 25) :: s7="<pseudo_partial_wave"
    character(len= 25) :: s9='<paw_radius'
    character(len= 25) :: s10='<blochl_local_io'
    character(len= 25) :: s11='<shape_function'
    character(len= 25) :: s12='<pseudo_valence_dens'
    inquire(file=p(it)%filename, exist=file_exists)
    if(.not.file_exists) then
        write(*,*) 'error : file not exist: ', p(it)%filename
        stop
    endif

    open(unit=98,file=p(it)%filename,action='read')
    write(*,*) 'now reading paw information for element type:',it
    write(*,*) 'from file: ',p(it)%filename

    call searchstring(s1)
    call get_atom_charge !get general information

    call searchstring(s9)
    call get_rcut !radius of augmentation sphere

    call searchstring(s2)
    call get_nstates_alloc !count #. states

    call searchstring(s4)
    call get_nr_alloc !count nr

    call searchstring(s2)
    call get_lstates_rcut !get l of each state and cutoff radius
    call get_nl
    call get_nm

    call searchstring(s4)
    call read_rr !read radial grid and derivatives
    call get_nrcut

    call read_p_phi !read projectors,ae and ps wf

    call searchstring(s10)
    call read_vloc  !read vlocal
    
    call searchstring(s11)
    call get_comp !read information on compensation charge
    call get_gl !generate shape fun
    call calc_qijlm !calculate qijlm

    call searchstring(s12)
    call read_ncoret !get ps core density 

    write(*,*) 'finish reading ',p(it)%filename
    write(*,*) '============================='
    write(*,*)
    close(98)
contains
    subroutine get_comp
        implicit none

        character(len=10) :: comp_type
        integer :: i,i1,i2
        real*8  :: w1,w2
        
        do i=1,len_trim(line)
           if (line(i:i) == '"') line(i:i) = " "
        end do

        whereis=index(line,'type')+len('type')
        read(line(1+whereis:),*) comp_type

        select case (trim(adjustl(comp_type)))
            case ('sinc')
                p(it)%flg_comp=1
            case ('gauss')
                p(it)%flg_comp=2
            case ('exp')
                p(it)%flg_comp=3
            case ('bessel')
                p(it)%flg_comp=4
            case default
                write(*,*) 'unknown shape function:', comp_type
                stop
        end select

        write(*,*) 'shape function type for compesation charge:',p(it)%flg_comp

        whereis=index(line,'rc')+len('rc')
        read(line(1+whereis:),*) p(it)%rcomp
        call interpolate(p(it)%rr,p(it)%nr,p(it)%rcomp,i1,i2,w1,w2)
        p(it)%nrcomp=i2
        write(*,*) 'radius parameter of compensation charge:',p(it)%rcomp
        write(*,*) 'the corresponding radial grid is: ',p(it)%nrcomp

        rewind 98
    end subroutine get_comp
    
    subroutine get_gl
        use m_paw_atom, only: atompaw_shapebes
        use m_paw_numeric, only : paw_jbessel
        implicit none

        integer :: i,nl,nrc,ir
        real*8  :: rc,arg
        real*8  :: norm,pi=3.14159265359d0
        real*8  :: bes_alpha(2),bes_q(2)
        real*8  :: jbes1,besp,bespp,jbes2

        nl=p(it)%nl
        rc=p(it)%rcomp
        nrc=p(it)%nrcomp

        allocate(p(it)%gl(nrc,p(it)%nl*2-1),stat=stat)
        if(stat/=0) stop 'gl alloc problem'
        norm=0d0

        do i=1,2*nl-1 
            select case (p(it)%flg_comp)
            case (1)
                if(i==1) then
                    p(it)%gl(1,i)=1d0
                    do ir=2,nrc
                        arg=pi*p(it)%rr(ir)/rc
                        p(it)%gl(ir,i)=(sin(arg)/arg)**2
                    enddo
                else
                    p(it)%gl(1,i)=0d0
                    do ir=2,nrc
                        arg=pi*p(it)%rr(ir)/rc
                        p(it)%gl(ir,i)=(sin(arg)/arg)**2*p(it)%rr(ir)**(i-1)
                    enddo
                endif
            case (2)
                stop 'gauss form not implemented yet!'
            case (3)
                stop 'exponential form not implemented yet!'
            case (4)
                call atompaw_shapebes(bes_alpha,bes_q,i-1,rc)
!                 write(*,*) 'il,alpha,q,rc',i,bes_alpha(1:2),&
!                     bes_q(1:2),rc
                 do ir=1,nrc
                    call paw_jbessel(jbes1,besp,bespp,i-1,&
                        0,bes_q(1)*p(it)%rr(ir))
                    call paw_jbessel(jbes2,besp,bespp,i-1,&
                        0,bes_q(2)*p(it)%rr(ir))
                    p(it)%gl(ir,i)=bes_alpha(1)*jbes1+bes_alpha(2)*jbes2
!                    write(1000+i,*) p(it)%rr(ir),p(it)%gl(ir,i)
                enddo
            end select

            norm=sum(p(it)%rr(1:nrc)**(i+1)*p(it)%dr(1:nrc)*&
                p(it)%gl(:,i))
            write(*,*) 'norm of gl:',i,norm
!            if(abs(norm-1d0)>1d-5 .and. p(it)%flg_comp==4) then
            if(p(it)%flg_comp==4) then
                write(*,*) 'Note:for bessel form norm should be 1'
            endif
            p(it)%gl(:,i)=p(it)%gl(:,i)/norm
        enddo
!stop

!        do i=1,nrc
!            write(800,*) p(it)%rr(i),p(it)%gl(i,:)
!        enddo
!        do i=1,2*nl-1
!            write(800+10*it+i,*) p(it)%gl(i,:) 
!            write(800+10*it+i,*) p(it)%gl(i,:) 
!        enddo
    end subroutine get_gl

    subroutine read_ncoret
        implicit none

        read(98,*) p(it)%ncoretilde
        p(it)%ncoretilde=p(it)%ncoretilde/sqrt(12.5663706144d0)
        rewind(98)

        write(*,*) 'int 4pi ncoretilde r^2 dr',&
            sum(p(it)%ncoretilde*p(it)%rr**2*p(it)%dr)*12.5663706144d0
    end subroutine read_ncoret

    subroutine read_vloc
        implicit none

        read(98,*) p(it)%vloc
        p(it)%vloc=p(it)%vloc*0.282094792
!        write(*,*) 'vloc',p(it)%vloc(1),p(it)%vloc(p(it)%nr)
        rewind(98)
    end subroutine read_vloc

    subroutine get_rcut
        implicit none

        integer :: i

        do i=1,len_trim(line)
           if (line(i:i) == '"') line(i:i) = " "
        end do
        whereis=index(line,'rc')+len('rc')
        read(line(1+whereis:),*) p(it)%rcut

        if(p(it)%rcut>xmax) then
            write(*,*) 'error:rcut>xmax'
            stop
        endif
        if(p(it)%rcut>ymax) then
            write(*,*) 'error:rcut>ymax'
            stop
        endif
        if(p(it)%rcut>zmax) then
            write(*,*) 'error:rcut>zmax'
            stop
        endif
        rewind(98)
    end subroutine get_rcut

    subroutine get_nrcut
        implicit none

        integer :: is,i1,i2,nrm !rr(ir)<r<rr(jr),jr=ir+1
        integer :: ns
        real*8  :: w1,w2

        ns=p(it)%nstates
    
        call interpolate(p(it)%rr,p(it)%nr,p(it)%rcut,i1,i2,w1,w2)
        p(it)%nrcut=i2
        call interpolate(p(it)%rr,p(it)%nr,p(it)%rcut+rpad_r,i1,i2,w1,w2)
        p(it)%nr_int=i2
        
        write(6,*) 'radius of augmentation sphere is: ',real(p(it)%rcut)
        write(6,*) 'the corresponding radial grid is: ',p(it)%nrcut
        write(6,*) 'nr_int: ',p(it)%nrcut
        write(6,*) ' dx, dy, dz ',real(dx),real(dy),real(dz)
        write(6,*) ' rpad ',real(rpad)

        p(it)%nrough(1)=ceiling((p(it)%rcut+rpad)/dx)+1
        p(it)%nrough(2)=ceiling((p(it)%rcut+rpad)/dy)+1
        p(it)%nrough(3)=ceiling((p(it)%rcut+rpad)/dz)+1
        p(it)%nrough=p(it)%nrough*2
        if(p(it)%nrough(1)>nx) then
            write(6,*)'it,nrough(1)>nx',it,p(it)%nrough(1),nx
            stop
        endif
        if(p(it)%nrough(2)>ny) then
            write(6,*)'it,nrough(2)>ny',it,p(it)%nrough(2),ny
            stop
        endif
        if(p(it)%nrough(3)>nz) then
            write(6,*)'it,nrough(3)>nz',it,p(it)%nrough(3),nz
            stop
        endif
        p(it)%nfine =p(it)%nrough*nfovnr
        write(*,*) 'nrough,nfine',p(it)%nrough,p(it)%nfine

        allocate(p(it)%bx(p(it)%nrough(1),p(it)%nfine(1)),&
            p(it)%by(p(it)%nrough(2),p(it)%nfine(2)),&
            p(it)%bz(p(it)%nrough(3),p(it)%nfine(3)),stat=stat)
        if(stat/=0) stop 'b alloc problem'

        call get_b
    end subroutine get_nrcut

    subroutine get_b
        implicit none

        real*8  :: xx1(p(it)%nrough(1)),yy1(p(it)%nrough(1)),y2a1(p(it)%nrough(1))
        real*8  :: xx2(p(it)%nrough(2)),yy2(p(it)%nrough(2)),y2a2(p(it)%nrough(2))
        real*8  :: xx3(p(it)%nrough(3)),yy3(p(it)%nrough(3)),y2a3(p(it)%nrough(3))
        integer :: i,j
        real*8  :: x,y

        p(it)%dxf(1)=dx*dble(p(it)%nrough(1)-1)/dble(p(it)%nfine(1)-1)

        do i=1,p(it)%nrough(1)
            x=(i-1)*dx
            xx1(i)=x
        enddo

        ! read under Eq. (C.5) for explanation
        ! using delta functions to generate the beta matrices
        do i=1,p(it)%nrough(1)
            yy1=0d0
            yy1(i)=1d0
            call spline_dn(xx1,yy1,p(it)%nrough(1),y2a1)
            do j=1,p(it)%nfine(1)
                x=(j-1)*p(it)%dxf(1)
                call splint_dn(xx1,yy1,y2a1,p(it)%nrough(1),x,y)
                p(it)%bx(i,j)=y*p(it)%dxf(1)/dx
            enddo
        enddo
        
!==========================        
        p(it)%dxf(2)=dy*dble(p(it)%nrough(2)-1)/dble(p(it)%nfine(2)-1)

        do i=1,p(it)%nrough(2)
            x=(i-1)*dy
            xx2(i)=x
        enddo

        do i=1,p(it)%nrough(2)
            yy2=0d0
            yy2(i)=1d0
            call spline_dn(xx2,yy2,p(it)%nrough(2),y2a2)
            do j=1,p(it)%nfine(2)
                x=(j-1)*p(it)%dxf(2)
                call splint_dn(xx2,yy2,y2a2,p(it)%nrough(2),x,y)
                p(it)%by(i,j)=y*p(it)%dxf(2)/dy
            enddo
        enddo

!==================        
        p(it)%dxf(3)=dz*dble(p(it)%nrough(3)-1)/dble(p(it)%nfine(3)-1)

        do i=1,p(it)%nrough(3)
            x=(i-1)*dz
            xx3(i)=x
        enddo

        do i=1,p(it)%nrough(3)
            yy3=0d0
            yy3(i)=1d0
            call spline_dn(xx3,yy3,p(it)%nrough(3),y2a3)
            do j=1,p(it)%nfine(3)
                x=(j-1)*p(it)%dxf(3)
                call splint_dn(xx3,yy3,y2a3,p(it)%nrough(3),x,y)
                p(it)%bz(i,j)=y*p(it)%dxf(3)/dz
            enddo
        enddo
        write(*,*) 'it,dxyz,dxyzf',it,dx,dy,dz,p(it)%dxf
        
    end subroutine get_b

    subroutine checkline(keyword,place)
        implicit none

        character(len=*) :: keyword,place

        read(98,'(A)') line
        whereis=index(line,keyword)
        if(whereis.eq.0) then
            write(*,*) 'error : wrong line ',place
            write(*,*) 'got: ',line
            stop
        endif

    end subroutine checkline

    subroutine get_nl
        implicit none

        integer :: i, min_val, max_val, ns, nl
        integer, allocatable :: unique(:),lstate(:)

        ns=p(it)%nstates

        allocate(unique(ns),lstate(ns),stat=stat)
        if(stat/=0) stop 'unique for counting l alloc problem '

        lstate=p(it)%lstate
        
        min_val = minval(lstate)-1
        max_val = maxval(lstate)

        i=0

        do while (min_val<max_val)
            i = i+1
            min_val = minval(lstate, mask=lstate>min_val)
            unique(i) = min_val
        enddo
        
        p(it)%nl=i
        nl=p(it)%nl

        write(*,*) 'number of distinct l: ', nl
        
        allocate(p(it)%lstate_diff(nl),stat=stat)
        if(stat/=0) stop 'lstate_diff alloc problem'

        p(it)%lstate_diff=unique(1:nl)
        
        write(*,*) "the distinct l's are: ", p(it)%lstate_diff

        deallocate(unique,lstate)
    end subroutine
    
    subroutine get_lstates_rcut
        implicit none

        integer :: i,j
        character(len=2) :: w1='l=',w2='f='

        do j=1,p(it)%nstates
            read(98,'(A)') line
            do i=1,len_trim(line)
               if (line(i:i) == '"') line(i:i) = " "
            end do

            whereis=index(line,trim(w1))+len(trim(w1))
            read(line(1+whereis:),*) p(it)%lstate(j)
            whereis=index(line,trim(w2))
            if(whereis>0) then
                whereis=whereis+3
                read(line(whereis:),*) p(it)%occ_l(j)
            else
                p(it)%occ_l(j)=0d0
            endif
        enddo

        write(*,*) 'l of valence states: ',p(it)%lstate
        write(*,*) 'occupation of valence states: ',p(it)%occ_l
    end subroutine get_lstates_rcut

    subroutine get_nm
        implicit none

        integer :: i,j,adv,ms,nmp

        p(it)%mstates=0
        do i=1,p(it)%nstates
            p(it)%mstates=p(it)%mstates+(2*p(it)%lstate(i)+1)
        enddo
        write(*,*) 'mstates: ',p(it)%mstates

        ms=p(it)%mstates

        allocate(p(it)%ms_ls(ms),p(it)%mstate(ms),stat=stat)
        if(stat/=0) stop 'mstates alloc problem'
!        allocate(p(it)%p3d(ms,nx,ny,nz),stat=stat)
!        if(stat/=0) stop 'p3d alloc problem'

        j=1
        do i=1,p(it)%nstates
            select case (p(it)%lstate(i))
            case (0)
                p(it)%ms_ls (j)=i
                p(it)%mstate(j)=0
                j=j+1
            case (1)
                p(it)%ms_ls (j:j+2)=i

                p(it)%mstate(j)  =-1
                p(it)%mstate(j+1)= 0
                p(it)%mstate(j+2)= 1
                j=j+3
                p(it)%occ_l(i)=p(it)%occ_l(i)/3d0
            case (2)
                p(it)%ms_ls(j:j+4)=i

                p(it)%mstate(j)  =-2
                p(it)%mstate(j+1)=-1
                p(it)%mstate(j+2)= 0
                p(it)%mstate(j+3)= 1
                p(it)%mstate(j+4)= 2
                j=j+5
                p(it)%occ_l(i)=p(it)%occ_l(i)/5d0
            case default
                stop 'l>2 not implemented yet'
            end select
        enddo
        
        write(*,*) 'm of valence states: ',p(it)%mstate
        write(*,*) 'correspondence between 3d and 1d states: ',p(it)%ms_ls
        write(*,*) 'adjusted occupation of states: ',p(it)%occ_l

    end subroutine get_nm

    subroutine read_p_phi
        implicit none

        integer :: i,j,l
        integer :: nr,ns
        
        nr=p(it)%nr_int
        ns=p(it)%nstates

        do i=1,ns
            call searchstring(s6)
            read(98,*) p(it)%phi(i,:)
            call searchstring(s7)
            read(98,*) p(it)%phitilde(i,:)
            call searchstring(s5)
            read(98,*) p(it)%ptilde(i,:)
!            call spline_dn(p(it)%rr,p(it)%ptilde(i,:),nr,p(it)%y2a_p(:,i))
        enddo
        
!        write(*,*) 'phi, for each state: '
!        do i=1,ns
!            write(*,*) i,p(it)%phi(i,1:3)
!        enddo
!
!        write(*,*) 'phitilde, for each state: '
!        do i=1,ns
!            write(*,*) i,p(it)%phitilde(i,1:3)
!        enddo
!
!        write(*,*) 'ptilde, for each state: '
!        do i=1,ns
!            write(*,*) i,p(it)%ptilde(i,1:3)
!        enddo

         do i=1,ns
            do j=1,ns
               p(it)%mat_sp(i,j)     =sum(p(it)%rr(1:nr)*p(it)%rr(1:nr)&
                   *p(it)%dr(1:nr)&  
                   *p(it)%ptilde(i,1:nr)*p(it)%ptilde(j,1:nr))

               p(it)%mat_t (i,j)     = &
                     sum(p(it)%rr(1:nr)*p(it)%rr(1:nr)*&
                     p(it)%dr(1:nr)*p(it)%phi(i,1:nr)&
                     *p(it)%phi(j,1:nr)) -&
                     sum(p(it)%rr(1:nr)*p(it)%rr(1:nr)*&
                     p(it)%dr(1:nr)*p(it)%phitilde(i,1:nr)&
                     *p(it)%phitilde(j,1:nr))
            enddo
        enddo
        call transform_pphi(it)

        do i=1,ns
!            call spline_dn(p(it)%rr,p(it)%ptilde1(i,:),nr,p(it)%y2a_p1(:,i))
        enddo
        rewind 98
    end subroutine read_p_phi

    subroutine read_rr
        implicit none

        call checkline('<values>','before reading rr') 

        read(98,*) p(it)%rr
        p(it)%rr=max(p(it)%rr,1d-8)

        write(*,*) 'max r for functions:',maxval(p(it)%rr)

        call checkline('</values>',    'after reading rr') 
        call checkline('<derivatives>','before reading dr') 

        read(98,*) p(it)%dr

        call checkline('</derivatives>','after reading dr')

    end subroutine read_rr

    subroutine get_nr_alloc
        implicit none

        character(len=10) :: w1='istart=', w2='iend='
        integer :: istart, iend
        integer :: i, nr, ns
        
        do i=1,len_trim(line)
           if (line(i:i) == '"') line(i:i) = " "
        end do

        whereis=index(line,trim(w1))+len(trim(w1))
        read(line(1+whereis:),*) istart

        whereis=index(line,trim(w2))+len(trim(w2))
        read(line(1+whereis:),*) iend

        p(it)%nr=iend-istart+1

        write(*,*) 'nr: ', p(it)%nr

        nr=p(it)%nr
        ns=p(it)%nstates 

        allocate(p(it)%rr(nr),p(it)%dr(nr),stat=stat)
        if(stat/=0) stop 'rr alloc problem'
        allocate(p(it)%y2a_p(nr,ns),stat=stat)
        if(stat/=0) stop 'y2a_p alloc problem'
        allocate(p(it)%y2a_p1(nr,ns),stat=stat)
        if(stat/=0) stop 'y2a_p1 alloc problem'
        allocate(p(it)%ptilde(ns,nr),stat=stat)
        if(stat/=0) stop 'ptilde alloc problem'
        allocate(p(it)%phitilde(ns,nr),stat=stat)
        if(stat/=0) stop 'phitilde alloc problem'
        allocate(p(it)%phi(ns,nr),stat=stat)
        if(stat/=0) stop 'phi alloc problem'
        allocate(p(it)%ptilde1(ns,nr),stat=stat)
        if(stat/=0) stop 'ptilde alloc problem'
        allocate(p(it)%phitilde1(ns,nr),stat=stat)
        if(stat/=0) stop 'phitilde alloc problem'
        allocate(p(it)%phi1(ns,nr),stat=stat)
        if(stat/=0) stop 'phi alloc problem'
        allocate(p(it)%vloc(nr),stat=stat)
        if(stat/=0) stop 'vloc alloc problem'
        allocate(p(it)%ncoretilde(nr),stat=stat)
        if(stat/=0) stop 'ncoretilde alloc problem'
        allocate(p(it)%mat_t(ns,ns),stat=stat)
        if(stat/=0) stop 'mat_t alloc problem'
        allocate(p(it)%mat_sp(ns,ns),stat=stat)
        if(stat/=0) stop 'mat_sp alloc problem'

        p(it)%ptilde=0d0
        p(it)%phi=0d0
        p(it)%phitilde=0d0
        p(it)%ptilde1=0d0
        p(it)%phi1=0d0
        p(it)%phitilde1=0d0
        p(it)%mat_t=0d0
        p(it)%mat_sp=0d0

        
        rewind(98)
    end subroutine get_nr_alloc

    subroutine get_nstates_alloc
        implicit none

        integer :: ns

        p(it)%nstates=0

        gdo : do while (.true.)
            
            read(98,'(A)',iostat=stat) line
            whereis=index(line,trim(s3))
            if(whereis.ne.0) then
                exit gdo
            elseif(stat<0) then
                write(*,*) 'error : string not found: ', s3
            endif

            p(it)%nstates=p(it)%nstates+1
        enddo gdo
       
        write(*,*) 'nstates: ',p(it)%nstates
        ns=p(it)%nstates 

        allocate(p(it)%lstate(ns),stat=stat)
        if(stat/=0) stop 'lstate alloc problem'
        allocate(p(it)%occ_l(ns),stat=stat)
        if(stat/=0) stop 'lstate alloc problem'

    end subroutine get_nstates_alloc

    subroutine searchstring(string)
        implicit none

        character(len=25) :: string
        
        gdo : do while (.true.)
            
            read(98,'(A)',iostat=stat) line
            whereis=index(line,trim(string))
            if(whereis.ne.0) then
                exit gdo
            elseif(stat<0) then
                write(*,*) 'error : string not found: ', string
            endif
        enddo gdo

    end subroutine searchstring

    subroutine get_atom_charge
        implicit none

        integer :: i
        character(len=10) :: w1="symbol=", w2="Z=", w3="core=", w4="valence="

        do i=1,len_trim(line)
           if (line(i:i) == '"') line(i:i) = " "
        end do

        whereis=index(line,trim(w1))+len(trim(w1))
        read(line(1+whereis:),*) p(it)%symbol

        whereis=index(line,trim(w2))+len(trim(w2))
        read(line(1+whereis:),*) p(it)%Zat

        whereis=index(line,trim(w3))+len(trim(w3))
        read(line(1+whereis:),*) p(it)%core

        whereis=index(line,trim(w4))+len(trim(w4))
        read(line(1+whereis:),*) p(it)%val
   

    end subroutine get_atom_charge

    subroutine calc_qijlm
        implicit none

        integer :: im,jm,il,jl,ll,mm,qmmax,ms,i,is,js
        real*8  :: nijl,cg,sq4pi=sqrt(12.5663706144d0)
        real*8, allocatable  :: integrand(:)

        allocate(integrand(p(it)%nr))

        qmmax=(2*p(it)%nl-1)**2
        ms=p(it)%mstates

        allocate(p(it)%qijlm(qmmax,ms,ms),stat=stat)
        if(stat/=0) stop 'qijlm alloc problem'
        allocate(p(it)%sij(ms,ms),stat=stat)
        if(stat/=0) stop 'sij alloc problem'
!        allocate(p(it)%ssqinvij(ms),stat=stat)
!        if(stat/=0) stop 'ssqinvij alloc problem'
        
        do i=1,1
            ll=floor(sqrt(dble(i-1)))
            mm=i-ll**2-ll-1
            do is=1,ms
                do js=1,ms
                    il=p(it)%ms_ls (is)
                    im=p(it)%mstate(is)
                    jl=p(it)%ms_ls (js)
                    jm=p(it)%mstate(js)

                    integrand=(p(it)%phi1(il,:)*p(it)%phi1(jl,:)- &
                        p(it)%phitilde1(il,:)*p(it)%phitilde1(jl,:))   &
                        *p(it)%rr**(ll+2)*p(it)%dr
                    integrand(1)=0d0
                    nijl=sum(integrand) !Eq (14)
!                    write(411,*) 'it,is,js,ll,nijl',it,is,js,ll,nijl
                    call cg3(ll,mm,p(it)%lstate(il),im,p(it)%lstate(jl),jm,cg)
!                    write(412,*) 'it,is,js,ll,mm,cg',it,is,js,ll,mm,cg
                    p(it)%qijlm(i,is,js)=nijl*cg!*sq4pi
!                    write(410,*) 'i,is,js,qijlm',i,is,js,p(it)%qijlm(i,is,js)
                enddo
            enddo
        enddo

        do is=1,ms
            do js=1,ms
                if (is/=js .and. abs(p(it)%qijlm(1,is,js))>1d-6) then
                    write(*,*) 'need to use paw files with diagonalized sij!!'
                    stop
                endif
            enddo
        enddo

        write(*,*) 'sij'
        do is=1,ms
            do js=1,ms
                p(it)%sij(is,js)=p(it)%qijlm(1,is,js)*sq4pi !becuase L=0
!                p(it)%sinvij(is)=1d0/(1d0+p(it)%sij(is))-1d0
!                p(it)%ssqinvij(is)=1d0/sqrt(1d0+p(it)%sij(is))-1d0 !note here
!                write(413,*) p(it)%sij(is,:)
            enddo
            if(p(it)%sij(is,is)<-1d0+toll) p(it)%sij(is,is)=-1d0+toll
            !overlap matrix sij, is 0th moment
            write(*,*) p(it)%sij(is,:)
        enddo

        do i=1,qmmax
            ll=floor(sqrt(dble(i-1)))
            mm=i-ll**2-ll-1
            do is=1,ms
                do js=1,ms
                    il=p(it)%ms_ls (is)
                    im=p(it)%mstate(is)
                    jl=p(it)%ms_ls (js)
                    jm=p(it)%mstate(js)

                    integrand=(p(it)%phi(il,:)*p(it)%phi(jl,:)- &
                        p(it)%phitilde(il,:)*p(it)%phitilde(jl,:))   &
                        *p(it)%rr**(ll+2)*p(it)%dr
                    integrand(1)=0d0
                    nijl=sum(integrand)
!                    write(411,*) 'it,is,js,ll,nijl',it,is,js,ll,nijl
                    call cg3(ll,mm,p(it)%lstate(il),im,p(it)%lstate(jl),jm,cg)
 !                   write(412,*) 'it,is,js,ll,mm,cg',it,is,js,ll,mm,cg
                    p(it)%qijlm(i,is,js)=nijl*cg!*sq4pi
                    write(410,*) 'i,is,js,qijlm',i,is,js,p(it)%qijlm(i,is,js)
                enddo
            enddo
        enddo
    end subroutine calc_qijlm

end subroutine readpaw 
