Module m_reload
#include 'Custom_Types.h'
    use GEODST_IO, only: GEODST_DATA
    use RTFLUX_IO, only: RTFLUX_DATA
    use NHFLUX_IO, only: NHFLUX_DATA
    use LABELS_IO, only: LABELS_DATA
    use PWDINT_IO, only: PWDINT_DATA
    use COMPXS_IO, only: COMPXS_DATA
    use NDXSRF_IO, only: NDXSRF_DATA
    use ZNATDN_IO, only: ZNATDN_DATA
    implicit none 
    integer, parameter :: max_num_fuel_regions = 5000
    integer, parameter :: max_num_prescribed   = 1000
    real(8) :: peak_tolerance   = 0.05D0     ! 5% relative tolerance for peak power limits
    real(8) :: eps_power = 1.0d-3  ! 0.1%, relative tolerance when comparing power/burnup values
    real(8), parameter :: eps_power_diff = 1.0d-4   ! 0.01%, absolute tolerance when comparing relative power change
    real(8), parameter :: eps_keff  = 1.0d-5  ! 1 pcm, absolute tolerance when comparing keff values
    real(8), parameter :: eps_dist  = 0.01d0  ! 1%, relative tolerance when comparing assembly distance
    real(8), parameter :: power_check_range = 0.1d0 ! 10%, candidate positions giving 0~10% higher estimated peak power than the minimum will be checked in DIF3D calculation when dif3d_check = True
    real(8), parameter :: burnup_check_range= 0.1d0 ! 10%, candidate positions having 0~10% lower discharge burnup than the maximum will be checked 
    integer, parameter :: num_objectives = 7
    Flag_Char,parameter:: defined_objectives(num_objectives) = &
    &   (/'assembly_power','peak_power','burnup','reactivity','power_shape','power_deep','reactivity_slow'/)
        ! range of objective_func (defined below):
        !   = 'assembly_power'  : peak assembly power will be minimized 
        !   = 'peak_power'      : peak power density will be minimized 
        !   = 'burnup'          : discharge burnup will be maximized
        !   = 'reactivity'      : added reactivity is maximized
        !   = 'power_shape'     : minimize relative deviation (infinite norm sense) from a reference power distribution
        !   = 'power_deep'      : reload at position with minimum relative assembly power reference condition
        !   = 'reactivity_slow' : (added reactivity - increased peaking factor) is maximized
    
    ! these tolerance parameter should be passed from reload_strategy.py later
    character(len=10), parameter :: dif3d_dir     = 'work_dif3d'
    character(len=12), parameter :: persent_dir   = 'work_persent'
    character(len=15), parameter :: reload_log    = 'reload.previous'
    character(len=13), parameter :: reload_manual = 'reload.manual'
    character(len=13), parameter :: ref_power     = 'reload.refpwr'

    ! &control input
    Flag_Char     :: flux_data        = 'RTFLUX'   ! 'RTFLUX' / 'NHFLUX'
    FileName_Char :: isotxs_file      = 'ISOTXS'
    FileName_Char :: burnup_data      = 'reload.burnup'
    FileName_Char :: rebus_output_file= 'rebus.out'
    Char8         :: fuel_regions(max_num_fuel_regions) = 'none'
    integer       :: num_reload_comp  = 1          ! number of fresh fuel compositions
    logical       :: null_refueling   = .false.    ! used to compute assembly worth change with depletion without refueling
    logical       :: gamma_heating    = .false.
    logical       :: fixed_peak_power = .true.     ! if False, whether to apply fixed limit of peak power
    logical       :: consider_keff    = .false.
    logical       :: consider_worth   = .true.
    logical       :: consider_burnup  = .true.
    logical       :: adjust_cycle_length = .true.   ! if True, will estimate cycle length for next depletion cycle after refueling
    logical       :: use_persent      = .false.     ! whether to run PERSENT for reacitivity perturbation calculation
    logical       :: dif3d_check      = .false.     ! whether to run DIF3D for final check of candidate reloading position
    logical       :: dif3d_check_all  = .false.     ! if True, force to check all candidates passed to reload_rundif3d, in effect only when dif3d_check is True
    logical       :: normalize_power  = .true.      ! whether to renormalize assembly powers after refueling, this is only used to choose reloading position
    logical       :: update_refpwr    = .false.     ! if True, the reference assembly powers are updated when the radial power peaking is reduced
    Flag_Char     :: objective_func   = 'assembly_power'    ! must be among defined_objectives, objective functions to be optimized to determine the final reloading position among candidates 
    Flag_Char     :: reload_refpwr    = 'bol'       ! 'bol', 'user', 'flat', default reference power distribution is obtained from BOL
    ! after screening out those violating hard constraints of power/burnup/reactivity etc.
    real(8)       :: peak_power_limit = -1.0       ! a fixed limit of peak power density, <0 means using the peak value at BOL as the limit value
    real(8)       :: assembly_power_limit = -1.0   ! a fixed limit of peak assembly power, <0 means using the peak value at BOL as the limit value
    real(4)       :: min_fast_energy  = 1.0D5      ! 100 keV, cutoff energy of fast neutron
    real(8)       :: min_burnup       = -1.0       ! minimum discharge burnup constraint
    ! this is used to narrow down candidate pool of reloading positions, used if objective_func /= 'burnup'
    ! default to -1.0 means that the average burnup of all fuel assemblies at EOC will be used as minimum discharge burnup
    real(8)       :: min_deltak       = -1.0       ! required minimum increase of keff after refueling
    ! this is used to narrow down candidate pool of reloading positions, used if objective_func /= 'reactivity'
    ! default to -1.0 means that the average delta_k at possible refueling positions will be used as min_deltak
    real(8)       :: max_deltak       = -1.0       ! allowed maximum increase of keff by refueling, used when consider_keff = True
    ! default to -1.0 means that max_deltak will be set to 0.002 (200 pcm) 
    real(8)       :: max_worth        = -1.0       ! allowed maximum worth of reloaded fresh fuel assembly
    ! default to -1.0 means that max_worth will be set to 0.002 (200 pcm)     
    real(8)       :: min_cycle_length = 1.0        ! allowed minimum cycle length in day, in effect when adjust_cycle_length is True
    real(8)       :: eps_cycle        = 0.1d0      ! default 10%, relative tolerance when comparing cycle length
    real(8)       :: target_eoc_keff  = 1.002      ! desired keff at eoc of each cycle, in effect when adjust_cycle_length is True
    real(8)       :: deprate_ratio    = 1.0        ! depletion rate correction factor, used to adjust the targeted eoc keff = boc_keff - deprate_ratio*(boc_keff - target_eoc_keff)
    integer       :: num_persent_case = 1          ! number of candidates for which reactivity worth is calculated in PERSENT
    ! for variable for 'reacvitiy_slow'
    real(8)       :: wgt_pekfac       = 1.0        ! weight of
    real(8)       :: conv_peak_rho    = 1.5e-3     ! conversion factor btn peaking factor and reactivity. Empirically, 1.5e-3 gives same selection with pekpwr/eocpwr - (pekfac - eocpekfac)
                                                   ! 1.5e-3 : 1 pcm reactivity = 0.66 percent peaking factor

    logical       :: fixed_sequence   = .FALSE.
    real(8)       :: margin_power_limit = -0.01    ! power limit margin for screening out, the candidates whose power > power_limit*(1-margin_power_limit) will be screened out
                                                   ! margin_power_limit is only used in screening out, not used actual peak power limit at the final phase
    ! used to give prescribed candidate range
    integer :: candidate_range_size = 0
    integer :: candidate_range(2,max_num_prescribed) = 0
    integer, pointer :: candidate_range_pos(:)  ! (candidate_range_size)

    ! default of 0 means to check all candidate positions after screening with other hard constraints 
    ! &material input
    type composition_data
        logical :: defined = .false.
        ISOTXS_Char :: name
        integer :: niso = 0
        !integer :: fiso
        ISOTXS_Char, pointer :: isonam(:)   ! isotope labels
        Real4, pointer :: density(:)  ! homogenized densities
    end type

    integer, parameter :: max_niso_per_comp = 500
    type (composition_data), pointer :: RFUEL(:)  ! (num_reload_comp or num_reload_comp+1 depending on consider_worth)
    ! if consider_worth = True, RFUEL(num_reload_comp+1) stores the composition for withdrawn assembly position

    ! measuring uncertainty effect of composition
    LOGICAL :: measure_uncertainty = .False.
    integer :: num_test_comp = 0
    type (composition_data), pointer :: TFUEL(:) ! fuel composition for test (power and reactivity)
    
    ! not input
    real(8) :: eps_keoc         ! absolute tolerance when comparing EOC keff values, determined by eps_cycle as eps_cycle * keff_loss_rate * min_cycle_length
    integer :: bot_active_node  ! assumed fuel nodes located at the same axial section in all fuel assemblies
    integer :: top_active_node  ! if not the case, need to record active node positions for individual assemblies
    integer :: num_xnode                    ! NINTI
    integer :: num_ynode                    ! NINTJ
    integer :: num_xynod                    ! num_xnode * num_ynode
    integer :: num_fuel_region              ! set in reload_input
    integer, pointer :: assembly_map(:,:)   ! (num_xnode, num_ynode), map of assembly types, 
    ! currently used to locate fuel assemblies, assembly outside problem domain is designated as type 0
    integer :: num_hex_sector  = 0          ! number of 60 degree sectors in the core domain (hexagonal geometry)
    integer, pointer :: hexagon_rp(:,:)     ! (2,num_xnode*num_ynode)   stores ring, position number for each assembly
    
    integer :: reload_type  = 1     ! type of reloaded assembly
    integer :: num_potential_site   ! number of possible locations for reloading, set in reload_potential_site
    integer, pointer :: fuel_pos(:) ! (num_potential_site) fuel assembly indices, set in reload_potential_site

    ! used to check distance among candidate refueling positions and previously reloaded assemblies
    logical :: disperse_reload      ! determined in reload_potential site by checking existence of reload_log file
    real(8) :: disperse_mindist     ! minimum distance (in unit of assembly pitch) of refueling position to recently reloaded assemblies
    integer :: num_previous_site    ! number of recently loaded assemblies to be checked (read from reload_log)
    integer, pointer :: previous_pos(:)   ! (num_previous_site) locations (assembly index) for previously loaded assemblies, earlier reloaded ordered before later reloaded 
    real(8), pointer :: disperse_dist(:)  ! (num_potential_site) minimum distance of refueling position to previously reloaded assemblies 

    ! used to handle manual setting of reloading position
    logical :: manual_reload
    integer :: manual_reload_pos(2)    ! location ((IRING, IPOS) for hex core or (IX,IY) for Cartesian core) for refueling, read from external file reload_manual
    !   later in reload_potential_site, manual_reload_pos is converted to (IX,IY) index in hex core case

    type kerma_data
        logical :: defined = .false.
        integer :: ncmp
        integer :: ngroup = 0
        integer :: ggroup = 0
        real(8), pointer :: nheat(:,:)  ! (ngroup, ncmp)
        real(8), pointer :: gheat(:,:)  ! (ggroup, ncmp)
    end type
    type burnup_data_type
        logical :: defined = .false.
        integer :: nreg
        Char8, pointer :: regnam(:)       ! (nreg) region names
        integer, pointer :: nburn(:)      ! (nreg) burned cycles
        real(8), pointer :: mass0(:)      ! (nreg) initial heavy metal mass
        real(8), pointer :: burnup(:)
    end type

    type (kerma_data) :: kerma
    real(4), allocatable :: grpstr(:)     ! energy group structure, get from ISOTXS
    real(8), allocatable :: node_vol(:)   ! (NINTK), node volume in an assembly, set in reload_power
    real(8), allocatable :: refpwr(:)     ! (num_xynod) reference assembly powers, only for fuel assemblies
    real(8), allocatable :: bocpwr(:)     ! (num_xynod) assembly powers from REBUS output, only for fuel assemblies 
    real(8), allocatable :: eocpwr(:)     ! (num_xynod) assembly powers before reloading
    real(8) :: eocpekpwr, eocpekdens      ! peak assembly power / power density before reloading
    real(8) :: eocpekfac                  ! assembly (radial) power peaking factor before reloading
    real(8) :: avgasspwr                  ! average assembly power of the core (total power / # of assemblies)
    real(8), allocatable :: bocden(:)     ! (num_xynod) peak power density in each assembly at BOC
    real(8), allocatable :: eocden(:)     ! (num_xynod) peak power density in each assembly at EOC
    integer, allocatable :: peknod(:)     ! (num_xynod) node position in each assembly showing peak power density
    real(8), allocatable :: asspwr(:,:)   ! (num_xynod, num_reload_comp) powers of reloaded assembly at different positions
    real(8), allocatable :: pekpwr(:,:)   ! (num_xynod, num_reload_comp) peak assembly powers with reloaded assembly at different positions
    real(8), allocatable :: pekfac(:,:)   ! (num_xynod, num_reload_comp) assembly (radial) peaking factor powers with reloaded assembly at different positions
    integer, allocatable :: pekass(:,:)   ! (num_xynod, num_reload_comp) position (GEODST ordering) of assembly having peak power with reloaded assembly
    real(8), allocatable :: pekden(:,:)   ! (num_xynod, num_reload_comp) peak (node-averaged) power density corresponding to different reloading positions
    integer, allocatable :: pekpos(:,:,:) ! (2,num_xynod, num_reload_comp) node position (IASS,IZ) showing peak power density corresponding to different reload positions

    type (burnup_data_type) :: region_burnup   ! used to store burnup data of each active region, obtained from rebus output
    integer, allocatable :: asscyc(:)     ! (num_xynod) assembly irradiation cycles
    real(8), allocatable :: assbrn(:)     ! (num_xynod) assembly burnup, MWD/MT
    real(8), allocatable :: assmas(:)     ! (num_xynod) assembly fuel mass, MT
    real(8), allocatable :: fluence(:,:,:) ! (NINTI, NINTJ, NINTK) node averaged fast fluence (n/cm^2)
    real(8) :: boc_peak_power_density     ! now is determined as the maximum node-averaged quantity
    real(8) :: boc_peak_assembly_power
    real(8) :: boc_keff
    real(8) :: eoc_keff
    real(8) :: cycle_length                ! cycle length for this cycle, read from ABURN file
    real(8), allocatable :: deltak(:,:)    ! (num_xynod, num_reload_comp+1) change in keff due to composition replacement of specific fuel assembly
    real(8), allocatable :: asswth(:,:)    ! (num_xynod, num_reload_comp) total assembly worth of reloaded fuel assembly
    real(8) :: mindk
    integer :: num_pert_cases              ! number of perturbations in one PERSENT run

    real(8) :: persent_time(3) ! timing PERSENT run (total, flux, integration)
    
    ! interface file data
    type (GEODST_DATA) :: GEODST
    type (RTFLUX_DATA) :: RTFLUX, boc_RTFLUX
    type (NHFLUX_DATA) :: NHFLUX
    type (LABELS_DATA) :: LABELS
    type (PWDINT_DATA) :: PWDINT, boc_PWDINT
    type (COMPXS_DATA) :: COMPXS    ! existing COMPXS at EOC
    !type (COMPXS_DATA) :: RFUELXS   ! macro XS of reloaded fuel
    type (NDXSRF_DATA) :: NDXSRF
    type (ZNATDN_DATA) :: ZNATDN

    save

End Module