from ic4popsyn import populations as pop
from ic4popsyn import tools

# Number of systems (10 systems are for backup)
Nbin = 2000
# create a population of binaries
binSanaMDS = pop.Binaries(Nbin, model='sana_eccM&DS', mass_ranges=[5.,150.], qmin=0.1, alphas=[-2.3]) # available: sana12 / sana_eccm&ds
SinglePop = pop.Binaries(Nbin, model='sana_eccM&DS', single_pop=True, mass_ranges=[5.,150.], alphas=[-2.3]) # available: sana12 / sana_eccm&ds


# save the population as input for SEVN
binSanaMDS.save_sevn_input('SEVNIC_BSE', 0.02, 0.02, 0.0, 0.0, 'end', 'zams', 'delayed', 'delayed', 'end')
binSanaMDS.save_sevn_input('SEVNIC_BSE_placeholder')
#SinglePop.save_sevn_input('SEVNIC_SSE', 0.02, 0.02, 0.0, 0.0, 'end', 'zams', 'delayed', 'delayed', 'end')
