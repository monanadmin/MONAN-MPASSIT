!> @file
!! @brief Geometria estruturada (lat-lon, projeções, e no futuro curvilínea).
!!
!! Passo 1.1 da refatoração: este módulo apenas "embrulha" um ESMF_Grid
!! já criado (hoje, por model_grid::define_target_grid). No passo 1.3 ele
!! ganha o construtor próprio (lat-lon por parâmetros), e na fase 3 a
!! versão curvilínea genérica (centros, cantos e máscara), para o MOM6.
!!
!! Regra da biblioteca de regrid: depende só de ESMF (e MPI). Não usa
!! model_grid, program_setup nem nenhum outro módulo do MPASSIT.
module structured_grid_mod

  use esmf

  implicit none

  private

  !> Grade estruturada decomposta entre os PETs.
  !! Os componentes são de uso interno da biblioteca (regrid_mod), exceto
  !! clb/cub, que o usuário consulta para dimensionar seus arrays locais.
  type, public :: structured_grid_t
     type(ESMF_Grid) :: grid                  !< grade ESMF
     integer         :: clb(2) = 0            !< limite inferior local (índices globais) nos centros
     integer         :: cub(2) = 0            !< limite superior local (índices globais) nos centros
     logical         :: owns_grid = .false.   !< .true. se este tipo criou a grade
  end type structured_grid_t

  public :: structured_grid_wrap

contains

  !> Embrulha um ESMF_Grid existente. Não copia nem assume a posse da grade.
  !! Supõe 1 DE por PET (como no MPASSIT).
  !!
  !! @param[out] this  geometria estruturada
  !! @param[in]  grid  grade ESMF já criada, com coordenadas no CENTER
  subroutine structured_grid_wrap(this, grid)

    type(structured_grid_t), intent(out) :: this
    type(ESMF_Grid),         intent(in)  :: grid

    integer :: rc

    call ESMF_GridGet(grid, localDE=0, staggerloc=ESMF_STAGGERLOC_CENTER, &
                      computationalLBound=this%clb, &
                      computationalUBound=this%cub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
                           line=__LINE__, file=__FILE__)) then
       print *, 'FATAL ERROR: structured_grid_wrap: ESMF_GridGet'
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    this%grid      = grid
    this%owns_grid = .false.

  end subroutine structured_grid_wrap

end module structured_grid_mod
