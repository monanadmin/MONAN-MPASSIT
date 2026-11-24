!> @file
!! @brief Geometria de origem: malha não estruturada do MPAS.
!!
!! Passo 1.1 da refatoração: este módulo apenas "embrulha" um ESMF_Mesh
!! já criado (hoje, por model_grid::define_input_grid). No passo 1.2 ele
!! ganha o construtor próprio que lê o arquivo de grade do MPAS.
!!
!! Regra da biblioteca de regrid: depende só de ESMF (e MPI). Não usa
!! model_grid, program_setup nem nenhum outro módulo do MPASSIT.
module mpas_mesh_mod

  use esmf

  implicit none

  private

  !> Malha MPAS decomposta entre os PETs.
  !! Os componentes são de uso interno da biblioteca (regrid_mod).
  type, public :: mpas_mesh_t
     type(ESMF_Mesh) :: mesh                  !< malha ESMF (células = elementos)
     integer         :: ncells_local = 0      !< nº de células deste PET
     logical         :: owns_mesh = .false.   !< .true. se este tipo criou a malha
  end type mpas_mesh_t

  public :: mpas_mesh_wrap

contains

  !> Embrulha um ESMF_Mesh existente. Não copia nem assume a posse da malha.
  !!
  !! @param[out] this  geometria MPAS
  !! @param[in]  mesh  malha ESMF já criada
  subroutine mpas_mesh_wrap(this, mesh)

    type(mpas_mesh_t), intent(out) :: this
    type(ESMF_Mesh),   intent(in)  :: mesh

    integer :: rc, n

    call ESMF_MeshGet(mesh, numOwnedElements=n, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                           line=__LINE__, file=__FILE__)) then
       print *, 'FATAL ERROR: mpas_mesh_wrap: ESMF_MeshGet'
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    this%mesh         = mesh
    this%ncells_local = n
    this%owns_mesh    = .false.

  end subroutine mpas_mesh_wrap

end module mpas_mesh_mod
