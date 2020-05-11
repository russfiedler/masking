Some tools to aid the choice of selection of the number of processors to be used and show how masking is done.

make sure there is a link to a file named ocean_grid.nc 

get_num_pes_masked_cice.f90 :Write out the processor distribution for a Round robin distribution in CICE
get_pe_number_mom.f90   : Write out the processor distribution for MOM5

These Ferret scripts hequire an ice history file with ice concentration

 pred_work_rr.jnl   : Estimate the amount of work done in the thermo step per processor for a masked round robin distribution.
                      This is assumed to be proportional to the ice coverage per processor.
                      A small number of PEs with excess work can result in a blowout in total time.

 work_demo.jnl      : Example driver script for the above with 3 examples for a layout with tile dimensioned 36x30

