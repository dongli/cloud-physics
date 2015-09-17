module constants

    implicit none

    real(8), parameter :: kb = 1.3806488d-23  ! Boltzmann's constant [J K-1].
    real(8), parameter :: NA = 6.022140857d23 ! Avogadro's number [mol-1].
    real(8), parameter :: Md = 28.9644d0 ! Dry air molecular weight [g mol-1].
    real(8), parameter :: Mw = 18.0160d0 ! Water molecular weight [g mol-1].
    real(8), parameter :: Rg = kb*NA ! Universal gas constant [J K-1 mol-1]
    real(8), parameter :: Rd = Rg/Md ! Dry air constant [J g-1 K-1].
    real(8), parameter :: Rv = Rg/Mw ! Water vapor constant [J g-1 K-1].
    real(8), parameter :: T_freeze = 273.15d0 ! Freezing point temperature of water [K].
    real(8), parameter :: T_triple = 273.16d0 ! Triple point temperature of water [K].
    real(8), parameter :: rhow = 1.0d0 ! Water density [g cm-3].

    ! "Constants" that is approximatively constant within an ordinary situation.
    real(8), parameter :: cpd = 1.005d0 ! Specific heat capacity of dry air at constant pressure [J g-1 K-1].
    real(8), parameter :: cvd = 0.718d0 ! Specific heat capacity of dry air at constant volume [J g-1 K-1].
    real(8), parameter :: cpv = 1.870d0 ! Specific heat capacity of water vapor at constant pressure [J g-1 K-1].
    real(8), parameter :: cvv = 1.410d0 ! Specific heat capacity of water vapor at constant volume [J g-1 K-1].

    real(8), parameter :: PI = atan(1.0d0)*4.0d0 ! Just pie.

    ! Unit conversion constants.
    real(8), parameter :: cal_to_J = 4.1868d0   ! [cal] -> [J].
    real(8), parameter :: hPa_to_J_cm_3 = 1.0d-4 ! [hPa] -> [J cm-3].
    real(8), parameter :: J_to_erg = 1.0d7 ! [J] -> [erg]

end module constants
