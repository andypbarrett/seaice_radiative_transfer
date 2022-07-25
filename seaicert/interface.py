"""Defines interface with ccsm3_sir_de.for fortran"""

"""
Output variables - taken from what is written to file object 6
        write(6,330) calday
        write(6,331) eccf
        write(6,332) rlat
        write(6,333) coszrs(i)
        write(6,334) loctim(i)
        write(6,335) sol  # Solar top of atmosphere insolation W/m2
        write(6,336) alb  # Solar TOA albedo
        write(6,337) sab  # Solar TOA absorbed W/m2
        write(6,338) fsnstm  # Solar absorbed atmosphere W/m2
        write(6,339) frs  # Solar absorbed surface W/m2
        write(6,340) clt  # Solar clear TOA absorbed W/m2
        write(6,341) clatm  # Solar clear atmosphere absorbed W/m2
        write(6,342) cls  # Solar clear surface absorbed W/m2
        write(6,343) scf  # Solar cloud forcing W/m2
        write(6,344) albc # Solar clear sky albedo
        write(6,345) flw  # Longwave net up TOA W/m2
        write(6,346) fla  # Longwave net surface up W/m2
        write(6,347) fld  # Longwave net down at surface W/m2
        write(6,348) clt  # Longwave clear outgoing W/m2
        write(6,349) cls  # Longwave clear net surface W/m2
        write(6,350) cfl  # Longwave cloud forcing W/m2
        write(6,351) cfn  # Net cloud forcing W/m2

        pmb is pmidm1/100
        cloud particle radii and ice fraction
        rel - r_liq
        rei - r_ice
        fice - ice fraction
          write(6,399) k,pmb,rel(i,k),rei(i,k),fice(i,k)

        Heating rates
        qrsday - qrs * 86400  Solar heating rate
        qrlday - qrl * 86400  Longwave cooling rate
          write(6,199) k,pmb,qrsday,qrlday,qrsday+qrlday

        write(6,380) sols(1)
        write(6,381) solsd(1)
        write(6,384) vsfdir
        write(6,386) (1.-vsfdir)
        write(6,382) soll(1)
        write(6,383) solld(1)
        write(6,385) nifdir
        write(6,387) (1.-nifdir)
        write(6,388) fsds
        write(6,361) vsfrac
        write(6,362) (1.-vsfrac)
        write(6,389) frs
        write(6,390) albsrf
        write(6,391) F_SW_vs
        write(6,392) F_SW_ni
        write(6,393) F_SW_srf
        write(6,394) I_vs
        write(6,395) I_ni  
        write(6,101) F_SW_vs*(1.-I_vs),F_SW_ni*(1.-I_ni),
            write(6,111) zd(k),Tri_vs(k),Tri_ni(k)
                write(6,322) k,Q_SW_vs,Q_SW_ni,
                write(6,325) k,Q_SW_vs,Q_SW_ni,
                  write(6,323) k,Q_SW_vs,Q_SW_ni,
                  write(6,325) k,Q_SW_vs,Q_SW_ni,
                  write(6,324) k,Q_SW_vs,Q_SW_ni,
                  write(6,325) k,Q_SW_vs,Q_SW_ni,
        write(6,111) zd(klevp),Tri_vs(klevp),Tri_ni(klevp)
        write(6,223) F_SW_ocn_vs,F_SW_ocn_ni,F_SW_ocn_vs+F_SW_ocn_ni
         write(6,*) ' day of year (1..365)  = ',dayyr(i)
         write(6,*) ' latitude (-90 to +90) = ',rlat(i)
            write(6,99) k   ,pmidm1(i,k),tm1(i,k),qm1(i,k),o3mmr(i,k)
         write(6,571) ps(i)
         write(6,572) co2mix
         write(6,573)    ts(i)
         write(6,574)    tg(i)
         write(6,575)    sndpth(i)
         write(6,576)    rhos(i)
         write(6,577)    rs(i)
         write(6,578)    hpnd(i)
         write(6,579)    R_pnd(i)
         write(6,580)    hice(i)
         write(6,581)    R_ice(i)

      write(6,*) ' cosine solar zenith angle    = ',mu0(1)

        write(6,591) albs(1)
        write(6,592) albsd(1)
        write(6,593) albl(1)
        write(6,594) albld(1)

              write(6,111) k,Fdirdn_vs(1,k),Fdirup_vs(1,k),
            write(6,112) k,Fdirdn_vs(1,k),Fdirup_vs(1,k),
        write(6,*) 
          write(6,112) k,Fdirdn_ni(1,k),Fdirup_ni(1,k),
      write(6,778)
             write(6,5) hs(i),hp(i)
        write(6,4321) ns,rupdir(1,0),rupdif(1,0)
c..                 write(6,1233) ns,exptdn(i,plevp)
c..                 write(6,1234) ns,(fluxdn(i,plevp)-exptdn(i,plevp))

"""

import ctypes

# These must be the same as in init.f
PLON = 1
PLEV = 18
NXPT = 1
plevp = PLEV + 1
plond = PLON + 1 + 2*NXPT


class OutputCom(ctypes.Structure):
    """Defines c-types for output common black used to return
       results"""
    _fields_ = [("asdir", ctypes.c_float * plond),
                ("asdif", ctypes.c_float * plond),
                ("aldir", ctypes.c_float * plond),
                ("aldif", ctypes.c_float * plond),
                ("F_SW_ocn_vs", ctypes.c_float),
                ("F_SW_ocn_ni", ctypes.c_float),
                ]


class InputCom(ctypes.Structure):
    """Defines c types for input common block that is
       used to set input parameters"""
    _fields_ = [("dayyr_in", ctypes.c_float),
                ("rlat_in", ctypes.c_float),
                ("lev_in", ctypes.c_int * PLEV),
                ("pmidm1_in", ctypes.c_float * PLEV),
                ("tm1_in", ctypes.c_float * PLEV),
                ("qm1_in", ctypes.c_float * PLEV),
                ("o3mmr_in", ctypes.c_float * PLEV),
                ("cld_in", ctypes.c_float * PLEV),
                ("clwp_in", ctypes.c_float * PLEV),
                ("ps_in", ctypes.c_float),
                ("co2mix_in", ctypes.c_float),
                ("ts_in", ctypes.c_float),
                ("tg_in", ctypes.c_float),
                ("sndpth_in", ctypes.c_float),
                ("rhos_in", ctypes.c_float),
                ("rs_in", ctypes.c_float),
                ("hpnd_in", ctypes.c_float),
                ("R_pnd_in", ctypes.c_float),
                ("hice_in", ctypes.c_float),
                ("R_ice_in", ctypes.c_float),
                ]


# Assign library and common blocks for interface
crmlib = ctypes.CDLL("../1D_dE_CCSM/libtest.so")  # need generic definition

# Common blocks
input_common = InputCom.in_dll(crmlib, "input_")
output_common = OutputCom.in_dll(crmlib, "output_")

# Alias for main program function
test = crmlib.test_parameters_

