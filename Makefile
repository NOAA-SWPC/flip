# flip Makefile

TOP = ..

SSRCS=CMINOR.f                       \
      CTIPE-int.f                    \
      ELECXS.f                       \
      FLIP_GRID.f                    \
      initialize_module_parameters.f \
      INIT-PROFILES.f                \
      KEMPRN.f                       \
      MINORA.f                       \
      module_tirosdata.f             \
      Neut_Heating.f                 \
      Photoel-Freqs.f                \
      Rates.f                        \
      RSDENA_EXB.f                   \
      RSJACA.f                       \
      RSLPSD.f                       \
      RSLPST.f                       \
      RSPE2B.f                       \
      RSPRIM.f                       \
      RSTEMC_EXB.f                   \
      sub_ionize_ipe.f               \
      sub_tiros_init.f               \
      sub_tiros_ionize.f

include $(TOP)/Makefile.common
