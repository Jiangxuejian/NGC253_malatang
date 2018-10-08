# Xuejian Jiang
#%% this script reads tables and join them, generate ratios and upper limits
# Sep 2018

from astropy.table import Table, join
import matplotlib.pyplot as plt
#import matplotlib.axes as ax
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame

dv1 = 10.2  # velocity resolution of CO32
dv2 = 19.8

#%% read files

origin='lower'
co32 = Table.read('co32_int_man.dat',format='ascii')
#co32.sort(['y','x'])
hcn0 = Table.read('hcn_int_man.dat',format='ascii')  # with several stare pixels.
#hcn0.sort(['y','x'])
hcn = hcn0[0:90]  # exclude four outer pixels
hcop0 = Table.read('hcop_int_man.dat',format='ascii')# with several stare pixels.
#hcop0.sort(['y','x'])
hcop = hcop0[0:90]  # exclude four outer pixels
co10_0 = Table.read('co10-45m.txt', format='ascii')
#co10_0.sort(['y','x'])
#---IR---
lir0 = Table.read('n253_psf_lir.dat',format='ascii') # luminosities, provided by Qinghua
fir0 = Table.read('n253_psf_flux.dat',format='ascii') # fluxes, provided by Qinghua
# off_x off_y l70 el70 l100 el100 l160 el160 lir elir pos_xpixel  pos_ypixel  flux_p70  err_flux_p70  flux_p100  err_flux_p100  flux_p160  err_flux_p100
lir0.rename_column('off_x','x')
lir0.rename_column('off_y','y')
fir0.rename_column('off_x','x')
fir0.rename_column('off_y','y')
ir = join(lir0,fir0, keys=['x','y'])


#%% -   HCN  --------------------
hcn['Int'].unit='K km/s'
hcn['rms'].unit='K'
#hcn_int = hcn['Int']/1000
hcn['Int']= hcn['Int']/1000
hcn['rms']= hcn['rms']/1000
for i in np.arange(0,len(hcn)):
    if hcn['W1'][i] == 'N/A' or hcn['W1'][i] == 'yes':
        hcn['W1'][i] =150
    if hcn['W2'][i] == 'N/A' or hcn['W2'][i] == 'yes':
        hcn['W2'][i] =350
hcn['line_win'] = np.float64(hcn['W2']) - np.float64(hcn['W1'])
hcn_total_win = 750
hcn['base_win'] = hcn_total_win - hcn['line_win']
hcn['err'] = hcn['rms'] * np.sqrt(hcn['line_win'] * dv1) * np.sqrt(1 + hcn['line_win']/hcn['base_win'])  #Gao+1996, Tan+2018
hcn['err'].unit= 'K km/s'

#%% CO32
co32['Int'] = co32['Int']/1000.
co32['rms'] = co32['rms']/1000.
co32['Int'].unit='K km/s'
co32['rms'].unit='K'
#%% HCO+  --------------------

hcop['Int'].unit='K km/s'
hcop['rms'].unit='K'
#hcop_int = hcop['Int']/1000.
hcop['Int'] = hcop['Int']/1000.
hcop['rms'] = hcop['rms']/1000.
hcop['line_win'] = np.float64(hcop['W2']) - np.float64(hcop['W1'])
hcop_total_win = 750
hcop['base_win'] = hcop_total_win - hcop['line_win']
hcop['err'] = hcop['rms'] * np.sqrt(hcop['line_win'] * dv1) * np.sqrt(1 + hcop['line_win']/hcop['base_win'])  #Gao+1996, Tan+2018
hcop['err'].unit= 'K km/s'

#%% join tables

co10_co32=join(co10_0, co32, keys=['x','y'], join_type='left',table_names=['CO10', 'CO32'],
            uniq_col_name='{table_name}_{col_name}')
hcn_hcop = join(hcn, hcop, keys=['x','y'], join_type='left',table_names=['HCN', 'HCOp'],
            uniq_col_name='{table_name}_{col_name}')
co10_co32_hcn_hcop = join(co10_co32, hcn_hcop, keys=['x','y'], join_type='outer')
longtable = join(co10_co32_hcn_hcop, ir, keys=['y','x'], join_type='left')
longtable.reverse()     # make x and y descending
longtable = Table(longtable, masked=True)  # convert to masked table
#%% - CO ------------------

co32_int = longtable['CO32_Int']
co32_line_win = longtable['CO32_W2'] - longtable['CO32_W1']
co32_total_win = 790
co32_base_win = co32_total_win - co32_line_win
co10_line_win = longtable['CO10_W2'] - longtable['CO10_W1']
co10_total_win = 560
co10_base_win = co10_total_win - co10_line_win
longtable['CO32_err'] = longtable['CO32_rms'] * np.sqrt(co32_line_win * dv1) * \
         np.sqrt(1 + co32_line_win/co32_base_win)   #Gao+1996, Tan+2018
longtable['CO10_err'] = longtable['CO10_rms'] * np.sqrt(co10_line_win * \
         longtable['dV']) * np.sqrt(1 + co10_line_win/co10_base_win)   #Gao+1996, Tan+2018

#%% - calculate separation --------------
# - central coordinate ------
c_pix_idx = 51 # index of the central pixel
c_pix = SkyCoord(longtable['ra'][c_pix_idx], longtable['dec'][c_pix_idx],unit=(u.radian, u.radian))
rotation=(51.07-90)*u.deg
c_pix.transform_to(SkyOffsetFrame(origin=c_pix,rotation=rotation))


#%% Sorai 2000  separation ~ 1.7''
#c2 = SkyCoord('00h45m5.63s', '-25d33m40.5s', frame='fk4', equinox='B1950')  # Sorai 2000
#c2_pix_idx = 5
#c2_pix = SkyCoord(co10['x'][c2_pix_idx], co10['y'][c2_pix_idx],unit=(u.arcsec, u.arcsec))
#co10['sep'] = co10['x']*u.arcsec #initialization of a new column
#co10['radec'] = SkyCoord(co10['x'], co10['y'],unit=(u.arcsec, u.arcsec))
#co10['radec_off'] = co10['radec'].transform_to(SkyOffsetFrame(origin=c2_pix,rotation=rotation))
#for i in np.arange(0,18):
#    c_tmp = SkyCoord(co10['x'][i], co10['y'][i],unit=(u.arcsec, u.arcsec))
#    if co10['radec_off'][i].lon > 0:
#        co10['sep'][i] = c_tmp.separation(c2_pix).arcsec
#    else:
#        co10['sep'][i] = c_tmp.separation(c2_pix).arcsec * -1


#%% relative coordinates --------------
longtable['sep'] = longtable['x']*u.arcsec #initialization of a new column
longtable['radec'] = SkyCoord(longtable['HCN_ra'], longtable['HCN_dec'],unit=(u.radian, u.radian))
longtable['radec'].transform_to(SkyOffsetFrame(origin=c_pix,rotation=rotation))
longtable['radec_off'] = longtable['radec'].transform_to(SkyOffsetFrame(origin=c_pix,rotation=rotation))

for i in np.arange(0,len(longtable)):
    if longtable['HCN_ra'].mask[i] == 1:
        longtable['HCN_ra'][i] = longtable['ra'][i]
        longtable['HCN_dec'][i] = longtable['dec'][i]
    c_tmp = SkyCoord(longtable['HCN_ra'][i], longtable['HCN_dec'][i],unit=(u.radian, u.radian))
    if longtable['radec_off'][i].lon > 0:
        longtable['sep'][i] = c_tmp.separation(c_pix).arcsec
    else:
        longtable['sep'][i] = c_tmp.separation(c_pix).arcsec * -1
#hcop['x']=hcn['x']
#hcop['y']=hcn['y']
# hcop['sep'] = hcn['sep']
# hcop['radec'] = hcn['radec']
# hcop['radec_off'] = hcn['radec_off']


longtable['r_pc'] = longtable['sep'] * 0.0122   # physical scale

#co['sep'] = hcn['sep']
#co10['sep'] = hcn['sep']
longtable['disk_flag'] = np.where(abs(longtable['radec_off'].lat.arcsec) < 11, 1, 0)
#co['disk_flag'] = longtable['disk_flag']
#hcop['disk_flag'] = co['disk_flag']


#%% - ratios --------------------
longtable['co32_to_co10'] = (longtable['CO32_Int'] / longtable['CO10_Int'])
longtable['co32_to_co10'].format = '5.5f'
longtable['co32_to_co10'].unit='--'
longtable['hcn_to_co32'] = (longtable['HCN_Int'] / longtable['CO32_Int'])
longtable['hcn_to_co32'].format = '5.5f'
longtable['hcn_to_co32'].unit='--'
longtable['hcn_to_co10'] = (longtable['HCN_Int'] / longtable['CO10_Int'])
longtable['hcn_to_co10'].format = '5.5f'
longtable['hcn_to_co10'].unit='--'
longtable['hcop_to_co32'] = (longtable['HCOp_Int'] / longtable['CO32_Int'])
longtable['hcop_to_co32'].format = '5.5f'
longtable['hcop_to_co32'].unit='--'
longtable['hcop_to_co10'] = (longtable['HCOp_Int'] / longtable['CO10_Int'])
longtable['hcop_to_co10'].format = '5.5f'
longtable['hcop_to_co10'].unit='--'
longtable['hcn_to_hcop'] = longtable['HCN_Int'] / longtable['HCOp_Int']
longtable['hcn_to_hcop'].format = '5.5f'
longtable['hcn_to_hcop'] = longtable['HCN_Int'] / longtable['HCOp_Int']
longtable['hcn_to_hcop'].format = '5.5f'

longtable['co32_to_co10_err'] = longtable['co32_to_co10'] * np.sqrt((longtable['CO32_err']/
         longtable['CO32_Int'])**2 + (longtable['CO10_err']/longtable['CO10_Int'])**2)
longtable['hcn_to_co10_err'] = longtable['hcn_to_co10'] * np.sqrt((longtable['HCN_err']/
         longtable['HCN_Int'])**2 + (longtable['CO10_err']/longtable['CO10_Int'])**2)
longtable['hcn_to_co32_err'] = longtable['hcn_to_co32'] * np.sqrt((longtable['HCN_err']/
         longtable['HCN_Int'])**2 + (longtable['CO32_err']/longtable['CO32_Int'])**2)
longtable['hcop_to_co10_err'] = longtable['hcop_to_co10'] * np.sqrt((longtable['HCOp_err']/
         longtable['HCOp_Int'])**2 + (longtable['CO10_err']/longtable['CO10_Int'])**2)
longtable['hcop_to_co32_err'] = longtable['hcop_to_co32'] * np.sqrt((longtable['HCOp_err']/
         longtable['HCOp_Int'])**2 + (longtable['CO32_err']/longtable['CO32_Int'])**2)

longtable['hcn_to_hcop_err'] = longtable['hcn_to_hcop'] * np.sqrt \
        ((longtable['HCN_err']/longtable['HCN_Int'])**2 + (longtable['HCOp_err']/longtable['HCOp_Int'])**2)
#co32_int_map=co32_int.reshape(7,13)
#co32_int_map=co32_int_map.transpose()

#%% HCN and HCO+ upper limits
# calculate 2 sigma upper limits:

# first deal with non-detections:
for i in np.arange(0,len(longtable)):
    if longtable['HCN_det'][i] == 0:
        longtable['hcn_to_co10'][i] = -99
        longtable['hcn_to_co32'][i] = -99
        longtable['hcn_to_hcop'][i] = -99
    if longtable['HCOp_det'][i] == 0:
        longtable['hcop_to_co10'][i] = -99
        longtable['hcop_to_co10'][i] = -99
        longtable['hcn_to_hcop'][i] = -99

#hcn_int_map=longtable['HCN_Int'].reshape(7,13)
##hcn_int_map=hcn_int_map.transpose()
#hcn_ratio = longtable['HCN_to_co']
#hcn_ratio_map = hcn_ratio.reshape(7,13)
#hcn_ratio_map = hcn_ratio_map.transpose()
hcnhcop_ratio = longtable['hcn_to_hcop']
#hcnhcop_ratio_map = hcnhcop_ratio.reshape(7,13)
#hcnhcop_ratio_map = hcnhcop_ratio_map.transpose()
hcn1 = longtable[np.where(longtable['HCN_det'] == 1)]  # positive detection
hcn2 = longtable[np.where(longtable['HCN_det'] == 2)]  # tentative detection
longtable['HCN_u_lim'] = longtable['HCN_rms'] * longtable['HCN_line_win'] * 2

# calculate 2 sigma upper limits:
#hcop_int_map=hcop_int.reshape(7,13)
#hcop_int_map=hcop_int_map.transpose()
#hcop_ratio = hcop['ratio_to_co']
#hcop_ratio_map = hcop_ratio.reshape(7,13)
#hcop_ratio_map = hcop_ratio_map.transpose()
hcop1 = longtable[np.where(longtable['HCOp_det'] == 1)]
hcop2 = longtable[np.where(longtable['HCOp_det'] == 2)]
longtable['HCOp_u_lim'] = longtable['HCOp_rms'] * longtable['HCOp_line_win'] * 2
#hcop2['u_lim_bar'] = 0.2
longtable['u_lim_bar'] = 0.2


# co32['lir']   =  hcn2['lir']
# co32['elir']  =  hcn2['elir']
# co32['f_70'] =  hcn2['flux_p70']
# co32['e_70'] =  hcn2['e_p70']
# co32['f_100'] =  hcn2['flux_p100']
# co32['e_100'] =  hcn2['e_p100']
# hcn['lir'] =  hcn2['lir']
# hcn['elir'] =  hcn2['elir']
# hcn['f_70']  =  hcn2['flux_p70']
# hcn['e_70']  =  hcn2['e_p70']
# hcn['f_100'] =  hcn2['flux_p100']
# hcn['e_100'] =  hcn2['e_p100']
# hcop['lir'] =  hcn2['lir']
# hcop['elir'] =  hcn2['elir']
# hcop['f_70']  =  hcn2['flux_p70']
# hcop['e_70']  =  hcn2['e_p70']
# hcop['f_100'] =  hcn2['flux_p100']
# hcop['e_100'] =  hcn2['e_p100']


#%% normalized intensity -------------------
longtable['CO10_norm_Int'] =  longtable['CO10_Int']/longtable['CO10_Int'].max()
longtable['CO32_norm_Int'] =  longtable['CO32_Int']/longtable['CO32_Int'].max()
longtable['HCN_norm_Int']= longtable['HCN_Int']/longtable['HCN_Int'].max()
longtable['HCOp_norm_Int']= longtable['HCOp_Int']/longtable['HCOp_Int'].max()
longtable['lir_norm']= longtable['lir']/longtable['lir'].max()


#%% Luminosity, SFR and Mass
fhcn = 354.223 # GHz
fhcop =356.447 # GHz
D_n253 = 3.5 # Mpc
z_n253 = 0.000811 # redshift
longtable['Lhcn']   = 3.25e7 * 24.4 * longtable['HCN_Int']  * fhcn**-2  * D_n253**2 * (1+z_n253)**-3
longtable['Lhcop'] = 3.25e7 * 24.4 * longtable['HCOp_Int'] * fhcop**-2 * D_n253**2 * (1+z_n253)**-3
longtable['Lhcn_err']   = 3.25e7 * 24.4 * longtable['HCN_err']  * fhcn**-2  * D_n253**2 * (1+z_n253)**-3
longtable['Lhcop_err'] = 3.25e7 * 24.4 * longtable['HCOp_err'] * fhcop**-2 * D_n253**2 * (1+z_n253)**-3
longtable['sfr'] = longtable['lir'] * 1.5e-10
longtable['sfr_err'] = longtable['elir'] * 1.5e-10


#%% conversion from L to M of hcn
bw70 =((
        60 * u.micron).to(u.Hz, equivalencies=u.spectral()) - (
        85 * u.micron).to(u.Hz, equivalencies=u.spectral())) # bandwidth of 70um in Hz 
bw100 =((
        85 * u.micron).to(u.Hz, equivalencies=u.spectral()) - (
        125 * u.micron).to(u.Hz, equivalencies=u.spectral())) # bandwidth of 100um in Hz 
I_FIR = longtable['flux_p70']*1e-23 * bw70/1e9 + \
        longtable['flux_p100']*1e-23 * bw100/1e9 /1e-17
G0 = 4 * np.pi * I_FIR / 1.6e-3 /2.7e-3  #The FUV field strength G0
longtable['Mhcn'] = longtable['Lhcn'] * 496* G0**-0.24      # conversion factor of M/L(HCN)
longtable['Mhcop'] = longtable['Lhcop'] * 689* G0**-0.24      # conversion factor of M/L(HCN)
longtable['Mhcn_err'] = longtable['Lhcn_err'] * 496* G0**-0.24      # conversion factor of M/L(HCN)
longtable['Mhcop_err'] = longtable['Lhcop_err'] * 689* G0**-0.24      # conversion factor of M/L(HCN)

#hcn['Mhcn'] = hcn['Lhcn'] * 10      # conversion factor of M/L(HCN)
#hcop['Mhcop'] = hcop['Lhcop'] * 10      # conversion factor of M/L(HCN)
longtable['sfe_hcn'] = longtable['sfr'] / longtable['Mhcn']
longtable['sfe_hcn_err'] = longtable['sfe_hcn'] * np.sqrt(
  (longtable['sfr_err']/longtable['sfr'])**2 + 
  (longtable['Mhcn_err']/longtable['Mhcn'])**2)
longtable['sfe_hcop'] = longtable['sfr'] / longtable['Mhcop']
longtable['sfe_hcop_err'] = longtable['sfe_hcop'] * np.sqrt(
  (longtable['sfr_err']/longtable['sfr'])**2 + 
  (longtable['Mhcop_err']/longtable['Mhcop'])**2)



longtable.write('table_all.dat', format='ascii', overwrite=True)

#hcn2['norm_Int'] = hcn2['Int']/hcn['Int'].max()
#hcop2['norm_Int'] = hcop2['Int']/hcop['Int'].max()
#
#co.write('co_new.dat', format='ascii', overwrite=True)
#co10.write('co10_new.dat', format='ascii', overwrite=True)
#hcn.write('hcn_new.dat', format='ascii', overwrite=True)
#hcop.write('hcop_new.dat', format='ascii', overwrite=True)
#hcn2.write('hcn_up_lim.dat', format='ascii', overwrite=True)
#hcop2.write('hcop_up_lim.dat', format='ascii', overwrite=True)

#%% --  figures --------------------
#  HCN / CO and HCO+/CO--------------------
#plt.close()
#fig=plt.figure()
#plt.scatter(hcop['sep'],  hcn['ratio_to_co'], edgecolors='red', c='none', label=r'  HCN (4-3)/CO (3-2)')
#plt.scatter(hcop['sep'], hcop['ratio_to_co'], edgecolors='blue',c='none', marker='^', label=r'HCO$^+$ (4-3)/CO (3-2)')
#plt.ylim(0,0.07)
#plt.xlim(45, -45)
#plt.legend()
#
#plt.xlabel(r'projected distance (arcsec)', size=16)
#plt.ylabel(r'$f_{\rm dense}$' , size=16)
#plt.title('Dense gas fraction profile in NGC 253')
#
#fig.savefig('ratio_profile.png')
#
#
#
##  HCO+ /  --------------------
#plt.close()
#fig=plt.figure()
#plt.scatter(hcop['sep'], hcn['hcn_to_hcop'], c='black')
#plt.ylim(0.1,3)
#plt.xlim(-1,25)
#
#plt.xlabel(r'projected distance (arcsec)', size=16)
#plt.ylabel(r'$R_{\rm HCN/HCO+}$' , size=16)
#plt.title(r'$R_{\rm HCN/HCO+}$ in NGC 253')
#
#fig.savefig('hcn_hcop_ratio_profile.png')
