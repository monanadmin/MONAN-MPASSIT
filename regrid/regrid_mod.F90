!> @file
!! @brief Núcleo da biblioteca de regrid: interpolação entre geometrias
!!        com interface de arrays Fortran.
!!
!! Uso típico:
!!
!!   use regrid_mod
!!   type(regridder_t) :: rg
!!   call rg%create(src, dst, method=REGRID_BILINEAR, missing=-9999.0_8)
!!   call rg%apply(t2m_in,  t2m_out)    ! 2D: (ncells_local)       -> (nx_local, ny_local)
!!   call rg%apply(temp_in, temp_out)   ! 3D: (ncells_local, nlev) -> (nx_local, ny_local, nlev)
!!   call rg%destroy()
!!
!! - Os pesos são calculados UMA vez, em create, e reaproveitados em todo
!!   apply, inclusive para campos 3D com qualquer número de níveis.
!! - Pontos de destino sem fonte (fora da malha de origem) recebem
!!   "missing". Não é preciso lista de pontos não mapeados: o destino é
!!   pré-preenchido com missing e o ESMF só zera e soma os pontos que
!!   recebem contribuição (zeroregion = ESMF_REGION_SELECT).
!! - Cada PET passa apenas a sua parte dos arrays:
!!     origem  : as células locais, na ordem da malha (mpas_mesh_t%ncells_local)
!!     destino : o bloco local dst%clb:dst%cub da grade estruturada
!!
!! Fase 1: origem = malha MPAS, destino = grade estruturada.
!! Validado só o método bilinear; conservativo e nearest entram na fase 2.
!!
!! Regra da biblioteca de regrid: depende só de ESMF (e MPI). Não usa
!! model_grid, program_setup nem nenhum outro módulo do MPASSIT.
module regrid_mod

  use esmf
  use mpas_mesh_mod,       only : mpas_mesh_t
  use structured_grid_mod, only : structured_grid_t

  implicit none

  private

  integer, parameter, public :: REGRID_BILINEAR = 1  !< interpolação bilinear
  integer, parameter, public :: REGRID_CONSERVE = 2  !< conservativo de 1ª ordem (fase 2)
  integer, parameter, public :: REGRID_NEAREST  = 3  !< vizinho mais próximo da origem (fase 2)

  !> Objeto de interpolação: guarda os pesos (route handle) entre uma
  !! origem e um destino fixos.
  type, public :: regridder_t
     private
     type(ESMF_RouteHandle) :: rh
     type(ESMF_Field)       :: src2d, dst2d
     type(ESMF_Mesh)        :: src_mesh
     type(ESMF_Grid)        :: dst_grid
     integer                :: ncells_src   = 0
     integer                :: dst_shape(2) = 0
     real(ESMF_KIND_R8)     :: missing      = -9999.0_ESMF_KIND_R8
     logical                :: ready        = .false.
   contains
     procedure          :: create   => regridder_create
     procedure          :: destroy  => regridder_destroy
     procedure, private :: apply_2d => regridder_apply_2d
     procedure, private :: apply_3d => regridder_apply_3d
     generic            :: apply    => apply_2d, apply_3d
  end type regridder_t

contains

  !> Calcula os pesos de interpolação da origem para o destino.
  !!
  !! @param[inout] this     objeto de interpolação
  !! @param[in]    src      malha de origem
  !! @param[in]    dst      grade de destino
  !! @param[in]    method   REGRID_BILINEAR, REGRID_CONSERVE ou REGRID_NEAREST
  !! @param[in]    missing  valor dos pontos sem fonte (padrão -9999.)
  subroutine regridder_create(this, src, dst, method, missing)

    class(regridder_t),      intent(inout)        :: this
    type(mpas_mesh_t),       intent(in)           :: src
    type(structured_grid_t), intent(in)           :: dst
    integer,                 intent(in)           :: method
    real(ESMF_KIND_R8),      intent(in), optional :: missing

    type(ESMF_RegridMethod_Flag) :: esmf_method
    integer                      :: rc, srcterm

    if (this%ready) call this%destroy()

    this%missing = -9999.0_ESMF_KIND_R8
    if (present(missing)) this%missing = missing

    select case (method)
    case (REGRID_BILINEAR)
       esmf_method = ESMF_REGRIDMETHOD_BILINEAR
    case (REGRID_CONSERVE)
       esmf_method = ESMF_REGRIDMETHOD_CONSERVE
    case (REGRID_NEAREST)
       esmf_method = ESMF_REGRIDMETHOD_NEAREST_STOD
    case default
       call regrid_abort('regridder_create: metodo de interpolacao desconhecido')
    end select

    this%src_mesh   = src%mesh
    this%dst_grid   = dst%grid
    this%ncells_src = src%ncells_local
    this%dst_shape  = dst%cub - dst%clb + 1

    ! Campos 2D-modelo: servem para calcular os pesos e para os applies 2D
    this%src2d = ESMF_FieldCreate(src%mesh, typekind=ESMF_TYPEKIND_R8, &
                                  meshloc=ESMF_MESHLOC_ELEMENT, &
                                  name='regrid_src2d', rc=rc)
    call check(rc, 'FieldCreate src2d')

    this%dst2d = ESMF_FieldCreate(dst%grid, typekind=ESMF_TYPEKIND_R8, &
                                  staggerloc=ESMF_STAGGERLOC_CENTER, &
                                  name='regrid_dst2d', rc=rc)
    call check(rc, 'FieldCreate dst2d')

    ! Mesmos parâmetros usados hoje no MPASSIT (srcTermProcessing=1),
    ! para que a saída seja idêntica à atual
    srcterm = 1
    call ESMF_FieldRegridStore(this%src2d, this%dst2d, &
                               regridmethod=esmf_method, &
                               routehandle=this%rh, &
                               srcTermProcessing=srcterm, &
                               unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                               rc=rc)
    call check(rc, 'FieldRegridStore')

    this%ready = .true.

  end subroutine regridder_create

  !> Interpola um campo 2D.
  !!
  !! @param[inout] this  objeto de interpolação (já criado)
  !! @param[in]    fin   valores nas células locais da origem (ncells_local)
  !! @param[out]   fout  valores no bloco local do destino (nx_local, ny_local)
  subroutine regridder_apply_2d(this, fin, fout)

    class(regridder_t), intent(inout) :: this
    real(ESMF_KIND_R8), intent(in)    :: fin(:)
    real(ESMF_KIND_R8), intent(out)   :: fout(:,:)

    real(ESMF_KIND_R8), pointer :: ps(:), pd(:,:)
    integer                     :: rc

    call require_ready(this, 'apply_2d')
    if (size(fin) /= this%ncells_src) &
       call regrid_abort('apply_2d: tamanho de fin difere do nº de celulas locais')
    if (any(shape(fout) /= this%dst_shape)) &
       call regrid_abort('apply_2d: forma de fout difere do bloco local do destino')

    call ESMF_FieldGet(this%src2d, farrayPtr=ps, rc=rc)
    call check(rc, 'FieldGet src2d')
    call ESMF_FieldGet(this%dst2d, farrayPtr=pd, rc=rc)
    call check(rc, 'FieldGet dst2d')

    ps = fin
    pd = this%missing

    call ESMF_FieldRegrid(this%src2d, this%dst2d, this%rh, &
                          zeroregion=ESMF_REGION_SELECT, rc=rc)
    call check(rc, 'FieldRegrid 2d')

    fout = pd

  end subroutine regridder_apply_2d

  !> Interpola um campo 3D, nível a nível, com os mesmos pesos do 2D.
  !! O ESMF permite aplicar um route handle a campos com dimensões não
  !! distribuídas (níveis) diferentes das usadas no cálculo dos pesos.
  !!
  !! @param[inout] this  objeto de interpolação (já criado)
  !! @param[in]    fin   valores na origem (ncells_local, nlev)
  !! @param[out]   fout  valores no destino (nx_local, ny_local, nlev)
  subroutine regridder_apply_3d(this, fin, fout)

    class(regridder_t), intent(inout) :: this
    real(ESMF_KIND_R8), intent(in)    :: fin(:,:)
    real(ESMF_KIND_R8), intent(out)   :: fout(:,:,:)

    type(ESMF_Field)            :: src3d, dst3d
    real(ESMF_KIND_R8), pointer :: ps(:,:), pd(:,:,:)
    integer                     :: rc, nlev

    call require_ready(this, 'apply_3d')
    nlev = size(fin, 2)
    if (size(fin, 1) /= this%ncells_src) &
       call regrid_abort('apply_3d: 1a dimensao de fin difere do nº de celulas locais')
    if (size(fout, 1) /= this%dst_shape(1) .or. size(fout, 2) /= this%dst_shape(2) &
        .or. size(fout, 3) /= nlev) &
       call regrid_abort('apply_3d: forma de fout incompativel com destino/niveis')

    src3d = ESMF_FieldCreate(this%src_mesh, typekind=ESMF_TYPEKIND_R8, &
                             meshloc=ESMF_MESHLOC_ELEMENT, &
                             ungriddedLBound=(/1/), ungriddedUBound=(/nlev/), &
                             name='regrid_src3d', rc=rc)
    call check(rc, 'FieldCreate src3d')
    dst3d = ESMF_FieldCreate(this%dst_grid, typekind=ESMF_TYPEKIND_R8, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             ungriddedLBound=(/1/), ungriddedUBound=(/nlev/), &
                             name='regrid_dst3d', rc=rc)
    call check(rc, 'FieldCreate dst3d')

    call ESMF_FieldGet(src3d, farrayPtr=ps, rc=rc)
    call check(rc, 'FieldGet src3d')
    call ESMF_FieldGet(dst3d, farrayPtr=pd, rc=rc)
    call check(rc, 'FieldGet dst3d')

    ps = fin
    pd = this%missing

    call ESMF_FieldRegrid(src3d, dst3d, this%rh, &
                          zeroregion=ESMF_REGION_SELECT, rc=rc)
    call check(rc, 'FieldRegrid 3d')

    fout = pd

    call ESMF_FieldDestroy(src3d, rc=rc)
    call ESMF_FieldDestroy(dst3d, rc=rc)

  end subroutine regridder_apply_3d

  !> Libera os pesos e os campos internos. Não destrói a malha nem a grade
  !! (elas pertencem a quem as criou).
  subroutine regridder_destroy(this)

    class(regridder_t), intent(inout) :: this
    integer :: rc

    if (.not. this%ready) return

    call ESMF_FieldRegridRelease(this%rh, rc=rc)
    call check(rc, 'FieldRegridRelease')
    call ESMF_FieldDestroy(this%src2d, rc=rc)
    call ESMF_FieldDestroy(this%dst2d, rc=rc)

    this%ready = .false.

  end subroutine regridder_destroy

  !-------------------------------------------------------------------------
  ! Utilitários internos
  !-------------------------------------------------------------------------

  subroutine require_ready(this, where)
    class(regridder_t), intent(in) :: this
    character(len=*),   intent(in) :: where
    if (.not. this%ready) call regrid_abort(where//': regridder nao foi criado (chame create)')
  end subroutine require_ready

  subroutine check(rc, msg)
    integer,          intent(in) :: rc
    character(len=*), intent(in) :: msg
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                           line=__LINE__, file=__FILE__)) &
       call regrid_abort('regrid_mod: '//msg)
  end subroutine check

  subroutine regrid_abort(msg)
    character(len=*), intent(in) :: msg
    print *, 'FATAL ERROR: ', trim(msg)
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  end subroutine regrid_abort

end module regrid_mod
