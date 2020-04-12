import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import colors
from scipy.interpolate import interp1d
import matplotlib.patches as patches
from scipy.signal import find_peaks
import json


from db import *
from mvmio import *

def add_timestamp(df, timecol='dt'):
  ''' add timestamp column assuming constant sampling in time '''
  df['timestamp'] = np.linspace( df.iloc[0,:][timecol] ,  df.iloc[-1,:][timecol] , len(df) )
  return df

def get_deltat(df, timestampcol='timestamp', timecol='dt'):
  ''' retrieve sampling time '''
  return df[timestampcol].iloc[2] - df[timestampcol].iloc[1]

def correct_sim_df(df):
  ''' Apply corrections to simulator data '''
  df['total_vol'] = df['total_vol'] * 0.1

def correct_mvm_df(df, pressure_offset=0, pv2_thr=50):
  ''' Apply corrections to MVM data '''
  # correct for miscalibration of pressure sensor if necessary
  df['pressure_pv1'] = df['pressure_pv1']  + pressure_offset
  df['airway_pressure'] = df['airway_pressure']  + pressure_offset

  #set MVM flux to 0 when out valve is open AND when the flux is negative
  df['flux'] = np.where ( df['out'] > pv2_thr , 0 , df['flux'] )
  df['flux'] = np.where ( df['flux'] < 0 , 0 , df['flux'] )

def apply_rough_shift(sim, mvm, manual_offset):
  imax1 = mvm[ ( mvm['dt']<8 ) & (mvm['dt']>2) ] ['flux'].idxmax()
  tmax1 = mvm['dt'].iloc[imax1]
  imax2 = sim[  (sim['dt']<8) & (sim['dt']>2)   ] ['total_flow'].idxmax()
  tmax2 = sim['dt'].iloc[imax2]
  shift  = tmax2 - tmax1
  if (manual_offset > 0 ) : print ("...Adding additional manual shift by [s]: ", manual_offset)
  mvm['dt'] = mvm['timestamp'] + shift  + manual_offset # add arbitrary shift to better match data

def apply_good_shift(sim, mvm, resp_rate, manual_offset):
  resp_period = 60./resp_rate

  sim_variable = 'total_flow' #'airway_pressure' #total_flow
  mvm_variable = 'flux'       #'airway_pressure' #flux
  sim['flux'] = np.where(sim[sim_variable]>0, sim[sim_variable], 0)
  sec = sim.dt[1]-sim.dt[0]
  first = sim.dt.to_list()[0]
  last  = sim.dt.to_list()[-1]
  sim_peaks = sim[(sim['dt']>first+5)&(sim['dt']<last-5)]#.rolling(1).mean()
  peak_rows, _ = find_peaks(sim_peaks['flux'].to_list(), prominence=sim_peaks['flux'].max()*0.5, distance=resp_period/sec*0.8)
  if (len(peak_rows)<2) :
    peak_rows, _ = find_peaks(sim_peaks['flux'].to_list(), prominence=sim_peaks['flux'].max()*0.3, distance=resp_period/sec*0.8)

  sim_peaks = sim_peaks.iloc[peak_rows]
  sim_peaks.sort_values(by=['dt'], inplace=True)
  sim_intervals = sim_peaks['dt'].diff().dropna().to_list()

  mvm['dt'] = mvm['timestamp']
  sec = mvm.dt[1]-mvm.dt[0]

  first = mvm.dt.to_list()[0]
  last  = mvm.dt.to_list()[-1]
  mvm_peaks = mvm[(mvm['dt']>first+5)&(mvm['dt']<last-5)]#.rolling(1).mean()
  peak_rows, _ = find_peaks(mvm_peaks[mvm_variable].to_list(), height=mvm_peaks[mvm_variable].max()*0.5, distance=resp_period/sec*0.8)
  if (len(peak_rows)<2) :
    peak_rows, _ = find_peaks(mvm_peaks[mvm_variable].to_list(), height=mvm_peaks[mvm_variable].max()*0.3, distance=resp_period/sec*0.8)
  mvm_peaks = mvm_peaks.iloc[peak_rows]
  mvm_peaks.sort_values(by=['dt'], inplace=True)
  mvm_intervals = mvm_peaks['dt'].diff().dropna().to_list()

  #print ("MVM PEAKS", mvm_peaks['dt'])
  #print ("SIM PEAKS", sim_peaks['dt'])

  mvm_peak_times = mvm_peaks['dt'].to_list()
  sim_peak_times = sim_peaks['dt'].to_list()
  mvm_peak_hgts = mvm_peaks[mvm_variable].to_list()
  sim_peak_hgts = sim_peaks['flux'].to_list()
  print ("I have identidied: ", len(mvm_peak_times), len(sim_peak_times))

  central_idx    = 30
  min_difference = 1e7
  tdiff          = 0
  for i in range (9) :
    offset = -4 + i
    subset_sim   = sim_peak_hgts[central_idx-10:central_idx+10]
    subset_mvm   = mvm_peak_hgts[central_idx-10 + offset :central_idx+10 + offset ]
    pdiff  = sum ( [ (x-y)*( x-y) for (x,y) in zip (subset_mvm, subset_sim)  ] )
    if pdiff < min_difference :
      min_difference = pdiff
      print ("minimisig distance: ",i, min_difference, mvm_peak_times[central_idx] - sim_peak_times [central_idx] )
      tdiff =  - np.mean( [ (x-y) for (x,y) in zip (mvm_peak_times, sim_peak_times)  ] )

  mvm_mean = np.nanmean(mvm_intervals)
  sim_mean = np.nanmean(sim_intervals)
  print(np.mean(mvm_intervals), np.mean(sim_intervals))

  interval = mvm_peaks.dt.to_list()[10]-sim_peaks.dt.to_list()[5]

  delay =  tdiff
  #delay = - ( interval - int(interval/mvm_mean)*mvm_mean )

  print('inspiratory rate ',60./mvm_mean, 'cycle / min')
  print('delay ',delay, 's')

  if (abs(manual_offset) > 0 ) : print (".. adding adidtional shift [s]", manual_offset )
  mvm['dt'] += delay + manual_offset

  """
  ax = mvm.plot( x='dt',y='flux')
  sim.plot(ax=ax, x='dt',y='flux')
  ax.plot(sim_peaks['dt'].to_list(),sim_peaks['flux'].to_list(),'x')
  ax.plot(mvm_peaks['dt'].to_list(),mvm_peaks['flux'].to_list(),'x')
  plt.show()
  """

def add_pv2_status(df):
  df['out_diff'] = df['out'].diff().fillna(0)
  df['out_status'] = np.where(df['out_diff'] == 0, 'steady', np.where(df['out_diff'] < 0, 'closing', 'opening'))

def get_start_times(df, thr=50, dist=0.1, quantity='out', timecol='dt'):
  ''' Get times where a given quantity starts being above threshold
  times_open = df[ ( df[quantity]>thr ) ][timecol]
  start_times = [ float(times_open.iloc[i]) for i in range(0, len(times_open)-1) if times_open.iloc[i+1]-times_open.iloc[i] > 0.1  or i == 0  ]
  '''
  start_times = df[df['out_status'] == 'closing']['dt'].unique()
  return start_times

def get_muscle_start_times(df, quantity='muscle_pressure', timecol='dt'):
  ''' True start time of breaths based on muscle_pressure (first negative time) '''
  negative_times = df[ ( df[quantity]<0 ) ][dt]    #array of negative times
  start_times = [ float (negative_times.iloc[i])  for i in range(0,len(negative_times)) if negative_times.iloc[i]-negative_times.iloc[i-1] > 1.  or i == 0  ]   #select only the first negative per bunch
  return start_times

def get_reaction_times(df, start_times, quantity='total_flow', timecol='dt', thr=5):
  ''' Reaction times based on total flow '''
  positive_flow_times = df[ ( df[quantity]>thr ) ][timecol]   # constant threshold to avoid artifact
  #flow_start_time = [ float (positive_flow_times.iloc[i])  for i in range(1,len(positive_flow_times)) if positive_flow_times.iloc[i]-positive_flow_times.iloc[i-1] > 1. ]   #select only the first negative per bunch
  flow_start_time = []

  for st in start_times :
    isFound = False
    for pf in positive_flow_times :
      if pf > st and isFound == False:
        isFound = True
        flow_start_time.append(pf)
    if isFound == False :
      flow_start_time.append(0)

  reaction_times = np.array ([100 * ( ft - st ) for ft,st in zip (flow_start_time,start_times) ])
  reaction_times = np.where ( abs(reaction_times)>1000, 0,reaction_times  )
  return reaction_times

def add_cycle_info(sim, mvm, start_times, reaction_times):
  ''' Add cycle start, reaction time, time-in-breath '''
  sim['start'] = 0
  sim['tbreath'] = 0
  sim['reaction_time'] = 0
  for s,f in zip (start_times,reaction_times) :
    sim.loc[sim.dt>s,'start']   = s
    sim.loc[sim.dt>s,'tbreath'] = sim.dt - s
    sim.loc[sim.dt>s,'reaction_time'] = f

  mvm['start'] = 0
  mvm['ncycle']= 0
  for i,s in enumerate(start_times) :
    mvm.loc[mvm.dt>s,  'start' ]   = s
    mvm.loc[mvm.dt>s,  'ncycle']   = i

def add_chunk_info(df):
  ''' Add mean values computed on simulator dataframe chunk '''
  cycle = df.groupby(['start']).mean()

  df['mean_pressure'] = 0

  for i,r in cycle.iterrows():
    df.loc[df.start == i, 'mean_pressure'] = r.total_flow

  df['norm_pressure']  = df['airway_pressure'] - df['mean_pressure']

  df['max_pressure'] = 0
  cycle = df.groupby(['start']).max()
  for i,r in cycle.iterrows():
    df.loc[df.start == i, 'max_pressure'] = r.total_flow

def add_clinical_values (df, max_R=250, max_C=100) :
  deltaT = get_deltat(df, timestampcol='dt')
  """Add for reference the measurement of "TRUE" clinical values as measured using the simulator"""

  #true resistance
  df['delta_pout_pin']        =  df['airway_pressure'] - df['chamber1_pressure']
  df['delta_vol']             = ( df['chamber1_vol'] * 2 ) .diff()
  df['airway_resistance']     =  df['delta_pout_pin'] / ( df['delta_vol'] / 60. )
  df.loc[ (abs(df.airway_resistance)>max_R) | (df.airway_resistance<0) ,"airway_resistance"] = 0

  #true compliance
  df['deltapin']     =  df['chamber1_pressure'].diff()
  df['delta_volin']  =  ( df['chamber1_vol']  + df['chamber2_vol']) . diff()
  df['compliance']   =  df['delta_volin']/df['deltapin']
  df.loc[abs(df.compliance)>max_C,"compliance"] = 0


def measure_clinical_values(df, start_times):
  ''' Compute tidal volume and other clinical quantities for MVM data '''
  deltaT = get_deltat(df)

  # determine inspiration end times
  inspiration_end_times   = df[df['out_status'] == 'opening']['dt'].unique()
  inspiration_start_times = start_times

  if (inspiration_end_times[0]< inspiration_start_times[0]) :
    inspiration_end_times = inspiration_end_times[1:]

# if df['out'].iloc[0] > 0: # if PV2 is open at beginning of considered dataset
#   tmp = [0]
#   for x in inspiration_end_times: tmp.append(x)
#   inspiration_end_times = tmp

  # propagate variables to all samples in each inspiration
  # i.e. assign to all samples in each cycle the same value of something
  df['is_inspiration'] = 0
  df['cycle_tidal_volume'] = 0
  df['cycle_peak_pressure'] = 0
  df['cycle_plateau_pressure'] = 0
  df['cycle_PEEP'] = 0
  df['cycle_Cdyn'] = 0
  df['cycle_Cstat'] = 0

  for i,(s,v) in enumerate(zip(inspiration_start_times, inspiration_end_times)):

    if i>=len(inspiration_start_times)-1 : continue

    this_inspiration   = (df.dt>s) & (df.dt<v) # all samples
    first_sample       = (df.dt == s) # start of inspiration
    last_sample        = (df.dt == v) # beginning of expiration
    next_inspiration_t = inspiration_start_times[i+1]

    df.loc[this_inspiration, 'is_inspiration'] = 1
    #measured inspiration
    df.loc[this_inspiration, 'cycle_tidal_volume']     = df[ this_inspiration ]['flux'].sum() * deltaT/60. * 100
    df.loc[this_inspiration, 'cycle_peak_pressure']    = df[ this_inspiration ]['airway_pressure'].max()
    df.loc[this_inspiration, 'cycle_plateau_pressure'] = df[ this_inspiration &( df.dt > v - 20e-3 ) & ( df.dt < v-10e-3 ) ]['airway_pressure'].mean()
    #not necessarily measured during inspiration
    df.loc[this_inspiration, 'cycle_PEEP']             = df[ ( df.dt > next_inspiration_t - 51e-3 ) & ( df.dt < next_inspiration_t+2e-3 ) ] ['airway_pressure'].mean()
    #print ("cycle_peak_pressure: " , df[ this_inspiration &( df.dt > v - 20e-3 ) & ( df.dt < v-10e-3 ) ]['airway_pressure'].mean() )

  # create cumulative variables (which get accumulated during a cycle)
  df['tidal_volume'] = 0
  for c in df['start'].unique():
    cycle_data = df[df['start']==c]
    cum = cycle_data['flux'].cumsum()
    df.loc[df.start == c, 'tidal_volume'] = cum

  # correct flow sum into volume (aka flow integral)
  df['tidal_volume'] *= deltaT/60.*100

  # set inspiration-only variables to zero outside the inspiratory phase
  df['tidal_volume'] *= df['is_inspiration']

  # calculate compliance and residence
  df['compliance'] = 0
  df['resistance'] = 0
  for c in df['start'].unique():
    df['vol_over_flux']           =  df['tidal_volume']*10. / df['flux'] * 60.
    df['delta_p_over_flux']       =  ( df['airway_pressure']-df['cycle_PEEP'])/ df['flux'] * 60.
    df['delta_vol_over_flux']     =  df['vol_over_flux'].diff()
    df['delta_delta_p_over_flux'] =  df['delta_p_over_flux'].diff()
    df['compliance']   =  df['delta_vol_over_flux']/ df['delta_delta_p_over_flux']*1000. # in ml/cmH2O
    df['resistance']   =  df['delta_p_over_flux']  - df['vol_over_flux']/df['compliance']*1000. # in

def add_run_info(df, dist=25):
  ''' Add run info based on max pressure '''
  df['mean_max_pressure'] = df['max_pressure'].rolling(4).mean()

  df['run'] = abs(df['mean_max_pressure'].diff())
  df['run'] = np.where(df['run'] > 2, 100, 0)

  starts  = df[df.run > 0]['dt'].values

  # add run number
  df['run'] = -1
  cc = 1
  for i in range(0,len(starts)-1):
    if starts[i+1] - starts[i] > dist:
      df.loc[(df['dt'] >= starts[i]) & (df['dt'] < starts[i+1]), 'run'] = cc
      cc += 1

  df['run'] = df['run']*10

def process_run(meta, objname, input_mvm, fullpath_rwa, fullpath_dta, columns_rwa, columns_dta, manual_offset=0., save=False, ignore_sim=False, mhracsv=None, pressure_offset=0, mvm_sep=' -> ', output_directory='plots_tmp'):
  # retrieve simulator data
  if not ignore_sim:
    df = get_simulator_df(fullpath_rwa, fullpath_dta, columns_rwa, columns_dta)
  else:
    print ("I am ignoring the simulator")

  # retrieve MVM data
  dfhd = get_mvm_df(fname=input_mvm, sep=mvm_sep)
  add_timestamp(dfhd)

  # apply corrections
  correct_mvm_df(dfhd, pressure_offset)
  correct_sim_df(df)

  #dfhd = dfhd[(dfhd['dt']>df['dt'].iloc[0]) & (dfhd['dt']<df['dt'].iloc[-1]) ]
  #rough shift for plotting purposes only - max in the first few seconds
  apply_rough_shift(sim=df, mvm=dfhd, manual_offset=manual_offset)
  #apply_good_shift(sim=df, mvm=dfhd, resp_rate=meta[objname]["Rate respiratio"], manual_offset=manual_offset)

  ##################################
  # cycles
  ##################################

  # add PV2 status info
  add_pv2_status(dfhd)

  # compute cycle start
  # start_times = get_muscle_start_times(df) # based on muscle pressure
  start_times    = get_start_times(dfhd) # based on PV2
  reaction_times = get_reaction_times(df, start_times)

  # add info
  add_cycle_info(sim=df, mvm=dfhd, start_times=start_times, reaction_times=reaction_times)

  ##################################
  # chunks
  ##################################

  add_chunk_info(df)

  # compute tidal volume etc
  add_clinical_values(df)
  measure_clinical_values(dfhd, start_times=start_times)

  #thispeep  = [ dfhd[dfhd.ncycle==i]['cycle_PEEP'].iloc[0] for
  measured_peeps      = []
  measured_volumes    = []
  measured_peak       = []
  measured_plateaus   = []
  real_tidal_volumes  = []
  real_plateaus       = []

  for i,nc in enumerate(dfhd['ncycle'].unique()) :

    this_cycle        = dfhd[ dfhd.ncycle==nc ]
    this_cycle_insp   = this_cycle[this_cycle.is_inspiration==1]

    if i > len(dfhd['ncycle'].unique()) -2 : continue
    #compute tidal volume in simulator df
    subdf             = df[ (df.dt>start_times[i]) & (df.dt<start_times[i+1]) ]

    subdf['total_vol_subtracted'] = subdf['total_vol'] - subdf['total_vol'].min()
    real_tidal_volume = subdf['total_vol_subtracted' ] .max()
    #compute plateau in simulator
    real_plateau = this_cycle_insp['airway_pressure'].iloc[-2]
    real_tidal_volumes.append(  real_tidal_volume   )
    real_plateaus.append (real_plateau)

    measured_peeps.append(  this_cycle['cycle_PEEP'].iloc[0])
    measured_volumes.append(this_cycle['cycle_tidal_volume'].iloc[0])
    measured_peak.append(   this_cycle['cycle_peak_pressure'].iloc[0])
    measured_plateaus.append(this_cycle['cycle_plateau_pressure'].iloc[0])

  """
  print ("measured_peeps", measured_peeps)
  print ("measured_volumes",measured_volumes)
  print ("measured_peak",measured_peak)
  print ("measured_plateau",measured_plateau)

  #plt.plot(dfhd['dt'], dfhd['out'], label='control out')
  plt.plot(dfhd['dt'], dfhd['airway_pressure'], label='pressure')
  plt.plot(dfhd['dt'], dfhd['is_inspiration'], label='is_inspiration')
  #plt.plot(dfhd['dt'], dfhd['cycle_tidal_volume'], label='cycle_tidal_volume')
  #plt.plot(dfhd['dt'], dfhd['cycle_PEEP'], label='cycle_PEEP')
  #plt.plot(dfhd['dt'], dfhd['tidal_volume'], label='tidal_volume')
  #plt.plot(dfhd['dt'], dfhd['cycle_tidal_volume'], label='cycle_tidal_volume')

  #plt.plot(df['dt'], df[''], label='')
  #plt.plot(df['dt'], df['tracheal_pressure'], label='true tracheal_pressure')
  #plt.plot(df['dt'], df['airway_pressure'], label='true airway_pressure')
  #plt.plot(df['dt'], df['chamber1_pressure'], label='true ch1_pressure')
  #plt.plot(df['dt'], df['chamber2_pressure'], label='true ch2_pressure')
  plt.plot(df['dt'], df['airway_resistance'], label='true airway_resistance')
  plt.plot(df['dt'], df['compliance'], label='true compliance')
  #plt.plot(dfhd['dt'], dfhd['compliance'], label='meas compliance')

  plt.legend()
  plt.show()
  print('Press ENTER to continue to real plots')
  input()
  """

  #raise RuntimeError('Valerio')

  ##################################
  # find runs
  ##################################

  add_run_info(df)

  ##################################
  # saving and plotting
  ##################################
  if save:
    df.to_hdf(f'{objname}.sim.h5', key='simulator')
    dfhd.to_hdf(f'{objname}.mvm.h5', key='MVM')

  if args.plot:
    colors = {  "muscle_pressure": "#009933"  , #green
      "sim_airway_pressure": "#cc3300" ,# red
      "total_flow":"#ffb84d" , #orange
      "tidal_volume":"#ddccff" , #orange
      "total_vol":"pink" , #
      "reaction_time" : "#999999", #
      "pressure" : "black" , #  blue
      "vent_airway_pressure": "#003399" ,# blue
      "flux" : "#3399ff" #light blue
    }

    #general service canavas
    ax = df.plot(x='dt', y='airway_pressure', label='airway_pressure [cmH2O]', c=colors['sim_airway_pressure'])
    df.plot(ax=ax, x='dt', y='total_flow',    label='total_flow      [l/min]', c=colors['total_flow'])
    df.plot(ax=ax, x='dt', y='run', label='run index')
    df.plot(ax=ax, x='dt', y='total_vol',          label='sim volume            [l/min]', c='red' )
    plt.plot(start_times, [0]*len(start_times), 'bo', label='real cycle start time')
    #df.plot(ax=ax , x='dt', y='muscle_pressure', label='muscle_pressure [cmH2O]', c=colors['muscle_pressure'], linewidth = linw)
    #df.plot(ax=ax, x='dt', y='breath_no', label='breath_no', marker='.')
    #df.plot(ax=ax, x='dt', y='tracheal_pressure', label='tracheal_pressure')
    #df.plot(ax=ax, x='dt', y='total_vol', label='total_vol')
    #plt.plot(df['dt'], df['total_vol']/10., label='total_vol [cl]')
    #dfhd.plot(ax=ax, x='dt', y='pressure', label='ventilator pressure [cmH2O]', c=colors['pressure'], linewidth = linw)
    dfhd.plot(ax=ax, x='dt', y='airway_pressure', label='ventilator airway pressure [cmH2O]', c=colors['vent_airway_pressure'])
    dfhd.plot(ax=ax, x='dt', y='flux',            label='ventilator flux            [l/min]', c=colors['flux'] )
    #dfhd.plot(ax=ax, x='dt', y='volume',          label='volume            [l/min]', c=colors['flux'] )
    #dfhd.plot(ax=ax, x='dt', y='out', label='out')
    #dfhd.plot(ax=ax, x='dt', y='in', label='in')
    #df.plot(ax=ax, x='dt', y='deriv_total_vol', label='deriv_total_vol [l/min]')
    #df.plot(ax=ax, x='dt', y='aux1', label='aux1')
    #df.plot(ax=ax, x='dt', y='aux2', label='aux2')
    #plt.gcf().suptitle(objname)
    #plt.plot(start_times, [0]*len(start_times), 'bo', label='real cycle start time')
    #df.plot(ax=ax, x='dt', y='reaction_time', label='reaction time [10 ms]', c=colors['reaction_time'])
    #plt.plot(   ,  100 *reaction_times,      label='reaction time ', marker='o', markersize=1, linewidth=0, c='red')
    #ax.
    for i,t in enumerate(start_times) :
      ax.text(t, 0.5, "%i"%i, verticalalignment='bottom', horizontalalignment='center', color='red', fontsize=14)

    ax.set_xlabel("Time [sec]")
    ax.legend(loc='upper center', ncol=2)

    #runs required for MHRA tests
    for i in range (len(meta)) :

      local_objname = "%s_%i"% ( objname[:-2] , i )

      PE = meta[local_objname]["Peep"]
      PI = meta[local_objname]["Pinspiratia"]
      RR = meta[local_objname]["Rate respiratio"]
      RT = meta[local_objname]["Resistance"]
      CM = meta[local_objname]["Compliance"]

      print ("Looking for R=%s, C=%s, RR=%s, PEEP=%s, PINSP=%s"%(RT,CM,RR,PE,PI) )

      my_selected_cycle = meta[local_objname]["cycle_index"] + 10

      print ("\nFor test [ %s ]  I am selecting cycle %i, starting at %f \n"%(meta[local_objname]["test_name"], my_selected_cycle , start_times[ my_selected_cycle ]))

      fig11,ax11 = plt.subplots()

      #make a subset dataframe for simulator
      dftmp = df[ (df['start'] >= start_times[ my_selected_cycle ] ) & ( df['start'] < start_times[ my_selected_cycle + 6])  ]
      #the (redundant) line below avoids the annoying warning
      dftmp = dftmp[ (dftmp['start'] >= start_times[ my_selected_cycle ] ) & ( dftmp['start'] < start_times[ my_selected_cycle + 6])  ]

      #make a subset dataframe for ventilator
      first_time_bin  = dftmp['dt'].iloc[0]
      last_time_bin   = dftmp['dt'].iloc[len(dftmp)-1]
      dfvent = dfhd[ (dfhd['dt']>first_time_bin) & (dfhd['dt']<last_time_bin) ]

      dftmp.loc[:, 'total_vol'] = dftmp['total_vol'] - dftmp['total_vol'].min()

      dftmp.plot(ax=ax11, x='dt', y='total_vol',         label='SIM tidal volume       [cl]', c=colors['total_vol'] , alpha=0.4)
      dftmp.plot(ax=ax11, x='dt', y='total_flow',        label='SIM flux            [l/min]', c=colors['total_flow'])
      dftmp.plot(ax=ax11, x='dt', y='airway_pressure',   label='SIM airway pressure [cmH2O]', c=colors['sim_airway_pressure'])

      dfvent.plot(ax=ax11,  x='dt', y='tidal_volume',    label='MVM tidal volume       [cl]', c=colors['tidal_volume'])
      dfvent.plot(ax=ax11,  x='dt', y='flux',            label='MVM flux            [l/min]', c=colors['flux'])
      dfvent.plot(ax=ax11,  x='dt', y='airway_pressure', label='MVM airway pressure [cmH2O]', c=colors['vent_airway_pressure'])

      ymin, ymax = ax11.get_ylim()
      ax11.set_ylim(ymin*1.45, ymax*1.55)
      ax11.legend(loc='upper center', ncol=2)
      title1="R = %i [cmH2O/l/s]         C = %2.1f [ml/cmH20]         PEEP = %s [cmH20]"%(RT,CM,PE )
      title2="Inspiration Pressure = %s [cmH20]       Frequency = %s [breath/min]"%(PI,RR)

      ax11.set_xlabel("Time [s]")

      xmin, xmax = ax11.get_xlim()
      ymin, ymax = ax11.get_ylim()
      ax11.text((xmax-xmin)/2.+xmin, 0.08*(ymax-ymin) + ymin,   title2, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      ax11.text((xmax-xmin)/2.+xmin, 0.026*(ymax-ymin) + ymin,  title1, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      nom_pressure = float(meta[local_objname]["Pinspiratia"])
      rect = patches.Rectangle((xmin,nom_pressure-2),xmax-xmin,4,edgecolor='None',facecolor='green', alpha=0.2)
      #add / remove nominal pressure band
      #ax11.add_patch(rect)

      nom_peep = float(meta[local_objname]["Peep"])
      rect = patches.Rectangle((xmin,nom_peep-0.1),xmax-xmin,0.5,edgecolor='None',facecolor='grey', alpha=0.3)
      #add / remove PEEP line
      #ax11.add_patch(rect)

      ax11.set_title ("Test n %s"%meta[objname]['test_name'], weight='heavy')
      figpath = "%s/%s_%s.pdf" % (output_directory, meta[objname]['Campaign'],  objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      print(f'Saving figure to {figpath}')
      fig11.savefig(figpath)

      dfvent.plot(ax=ax11,  x='dt', y='pressure_pv1', label='MVM PV1 pressure [cmH2O]', c='black')
      figpath = "%s/%s_service_%s.pdf" % (output_directory, meta[objname]['Campaign'],  objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      ax11.legend(loc='upper center', ncol=2)
      fig11.savefig(figpath)

      ####################################################
      '''try 30 cycles'''
      ####################################################
      """
      fig30c,ax30cycles = plt.subplots(3,1)
      ax30cycles = ax30cycles.flatten()

      my_selected_cycle = 10

      print (len(start_times) )
      for ii in range (3) :
        if len (start_times) < (ii+1) *10 : continue
        print ( my_selected_cycle + ii*10)
        #make a subset dataframe for simulator
        print (    start_times[ my_selected_cycle + ii*10 ] , start_times[ my_selected_cycle + ii*10 +11 ] )
        dftmp = df[ (df['start'] >= start_times[ my_selected_cycle + ii*10 ] ) & ( df['start'] < start_times[ my_selected_cycle + 11 + ii*10 ])  ]
        #the (redundant) line below avoids the annoying warning
        dftmp = dftmp[ ( dftmp['start'] >= start_times[ my_selected_cycle + ii*10  ] ) & ( dftmp['start'] < start_times[ my_selected_cycle + 11 + ii*10 ])  ]
        print (len (dftmp) , dftmp, dftmp['dt'].iloc[len(dftmp)-1] )
        #make a subset dataframe for ventilator
        first_time_bin  = dftmp['dt'].iloc[0]
        last_time_bin   = dftmp['dt'].iloc[len(dftmp)-1]
        dfvent = dfhd[ (dfhd['dt']>first_time_bin) & (dfhd['dt']<last_time_bin) ]

        dftmp.loc[:, 'total_vol'] = dftmp['total_vol'] - dftmp['total_vol'].min()
        dftmp.plot(ax=ax30cycles[ii], x='dt', y='total_flow',        label='SIM flux            [l/min]', c=colors['total_flow'] ,  legend=False)
        dftmp.plot(ax=ax30cycles[ii], x='dt', y='airway_pressure',   label='SIM airway pressure [cmH2O]', c=colors['sim_airway_pressure'] ,  legend=False)
        dfvent.plot(ax=ax30cycles[ii],  x='dt', y='flux',            label='MVM flux            [l/min]', c=colors['flux']  ,  legend=False)
        dfvent.plot(ax=ax30cycles[ii],  x='dt', y='airway_pressure', label='MVM airway pressure [cmH2O]', c=colors['vent_airway_pressure']   ,  legend=False )
        ax30cycles[ii].set_xlabel("Time [s]")


      #ymin, ymax = ax30cycles.get_ylim()
      #ax30cycles.set_ylim(ymin*1.4, ymax*1.5)
      #ax30cycles.legend(loc='upper center', ncol=2)
      #title1="R = %i [cmH2O/l/s]         C = %i [ml/cmH20]         PEEP = %s [cmH20]"%(RT,CM,PE )
      #title2="Inspiration Pressure = %s [cmH20]       Frequency = %s [breath/min]"%(PI,RR)

      '''
      xmin, xmax = ax30cycles.get_xlim()
      ymin, ymax = ax30cycles.get_ylim()
      ax30cycles.text((xmax-xmin)/2.+xmin, 0.08*(ymax-ymin) + ymin,   title2, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      ax30cycles.text((xmax-xmin)/2.+xmin, 0.026*(ymax-ymin) + ymin,  title1, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      nom_pressure = float(meta[local_objname]["Pinspiratia"])
      rect = patches.Rectangle((xmin,nom_pressure-2),xmax-xmin,4,edgecolor='None',facecolor='green', alpha=0.2)
      ax30cycles.add_patch(rect)

      nom_peep = float(meta[local_objname]["Peep"])
      rect = patches.Rectangle((xmin,nom_peep-0.1),xmax-xmin,0.5,edgecolor='None',facecolor='grey', alpha=0.3)
      ax30cycles.add_patch(rect)
      '''
      #ax30cycles.set_title ("Test n %s"%meta[objname]['test_name'])
      figpath = "%s/%s_30cycles_%s.pdf" % (output_directory, meta[objname]['Campaign'],  objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      print(f'Saving figure to {figpath}')
      fig30c.savefig(figpath)
      """

      ####################################################
      '''plot the avg wfs'''
      ####################################################

      fig2, ax2 = plt.subplots()
      #make a subset dataframe for simulator
      dftmp = df[ (df['start'] >= start_times[ 4 ] ) & ( df['start'] < start_times[ min ([35,len(start_times)-1] )  ])]
      dftmp['dtc'] = df['dt'] - df['start']

      #make a subset dataframe for ventilator
      first_time_bin  = dftmp['dt'].iloc[0]
      last_time_bin   = dftmp['dt'].iloc[len(dftmp)-1]
      dfvent = dfhd[ (dfhd['dt']>first_time_bin) & (dfhd['dt']<last_time_bin) ]
      dfvent['dtc'] = dfvent['dt'] - dfvent['start']
      dfvent = dfvent.sort_values('dtc')

      dftmp.loc[:, 'total_vol'] = dftmp['total_vol'] - dftmp['total_vol'].min()

      dftmp.plot(ax=ax2, x='dtc', y='total_vol',         label='SIM tidal volume       [cl]', c=colors['total_vol'] ,          marker='o', markersize=0.3, linewidth=0)
      dftmp.plot(ax=ax2, x='dtc', y='total_flow',        label='SIM flux            [l/min]', c=colors['total_flow'],          marker='o', markersize=0.3, linewidth=0)
      dftmp.plot(ax=ax2, x='dtc', y='airway_pressure',   label='SIM airway pressure [cmH2O]', c=colors['sim_airway_pressure'], marker='o', markersize=0.3, linewidth=0)

      dfvent.plot(ax=ax2,  x='dtc', y='tidal_volume',    label='MVM tidal volume       [cl]', c=colors['tidal_volume'],         marker='o', markersize=0.3, linewidth=0.2)
      dfvent.plot(ax=ax2,  x='dtc', y='flux',            label='MVM flux            [l/min]', c=colors['flux'],                 marker='o', markersize=0.3, linewidth=0.2)
      dfvent.plot(ax=ax2,  x='dtc', y='airway_pressure', label='MVM airway pressure [cmH2O]', c=colors['vent_airway_pressure'], marker='o', markersize=0.3, linewidth=0.2)

      ymin, ymax = ax2.get_ylim()
      ax2.set_ylim(ymin*1.4, ymax*1.5)
      ax2.legend(loc='upper center', ncol=2)
      title1="R = %i [cmH2O/l/s]         C = %2.1f [ml/cmH20]        PEEP = %s [cmH20]"%(RT,CM,PE )
      title2="Inspiration Pressure = %s [cmH20]       Frequency = %s [breath/min]"%(PI,RR)

      ax2.set_xlabel("Time [s]")

      xmin, xmax = ax2.get_xlim()
      ymin, ymax = ax2.get_ylim()
      ax2.text((xmax-xmin)/2.+xmin, 0.08*(ymax-ymin) + ymin,   title2, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      ax2.text((xmax-xmin)/2.+xmin, 0.026*(ymax-ymin) + ymin,  title1, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      nom_pressure = float(meta[local_objname]["Pinspiratia"])
      rect = patches.Rectangle((xmin,nom_pressure-2),xmax-xmin,4,edgecolor='None',facecolor='green', alpha=0.2)
      ax2.add_patch(rect)

      nom_peep = float(meta[local_objname]["Peep"])
      rect = patches.Rectangle((xmin,nom_peep-0.1),xmax-xmin,0.5,edgecolor='None',facecolor='grey', alpha=0.3)
      ax2.add_patch(rect)

      ax2.set_title ("Test n %s"%meta[objname]['test_name'])
      figpath = "%s/%s_avg_%s.png" % (output_directory, meta[objname]['Campaign'] , objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      print(f'Saving figure to {figpath}')
      fig2.savefig(figpath)

      ####################################################
      '''summary plots'''
      ####################################################

      #correct for outliers ?
      measured_peeps      = measured_peeps[3:-3]
      measured_plateaus   = measured_plateaus[3:-3]
      measured_peak       = measured_peak[3:-3]
      measured_volumes    = measured_volumes[3:-3]

      mean_peep    = np.mean(measured_peeps)
      mean_plateau = np.mean(measured_plateaus)
      mean_peak    = np.mean(measured_peak)
      mean_volume  = np.mean(measured_volumes)
      rms_peep    = np.std(measured_peeps)
      rms_plateau = np.std(measured_plateaus)
      rms_peak    = np.std(measured_peak)
      rms_volume  = np.std(measured_volumes)

      figs,axes = plt.subplots(2,2)
      axs = axes.flatten()
      #axs.set_title("PEEP", "", "a", "")
      nom_peep_low = nom_peep - 2 - 0.04 * nom_peep
      nom_peep_wid = 4 + 0.08 * nom_peep
      axs[0].hist ( measured_peeps  , bins=50,  range=(  min([ mean_peep,nom_peep] )*0.6 , max( [mean_peep,nom_peep] ) *1.4  )   )
      aa = patches.Rectangle( (nom_peep_low, axs[0].get_ylim()[0]  ) , nom_peep_wid , axs[0].get_ylim()[1] , edgecolor='red' , facecolor='green' , alpha=0.2)
      axs[0].add_patch(aa)
      axs[0].set_title("PEEP [cmH20], nominal: %i [cmH20]"%nom_peep)

      nominal_plateau = meta[objname]["Pinspiratia"]
      nominal_plateau_low = nominal_plateau - 2 - 0.04 * nominal_plateau
      nominal_plateau_wid = 4 + 0.08 * nominal_plateau
      axs[1].hist ( measured_plateaus, bins=100, range=(   min([ mean_plateau,nominal_plateau] )*0.8 , max( [mean_plateau,nominal_plateau] ) *1.3  )   )
      aa = patches.Rectangle( (nominal_plateau_low, axs[0].get_ylim()[0]  ) , nominal_plateau_wid , axs[0].get_ylim()[1] , edgecolor='red' , facecolor='green' , alpha=0.2)
      axs[1].add_patch(aa)
      axs[1].set_title("plateau [cmH20], nominal: %s [cmH20]"%nominal_plateau)

      #peak - put plateau instead
      #axs[2].hist ( measured_peak   , bins=100, range=(mean_peak*0.8 , mean_peak*1.3)  )
      #axs[2].set_title("peak pressure [cmH20]")

      nominal_plateau     =   np.nanmean( real_plateaus)   # float ( meta[objname]["Tidal plateau"] ) / 10.
      nominal_plateau_low = nominal_plateau - 2 - 0.04 * nominal_plateau
      nominal_plateau_wid = 4 + 0.08 * nominal_plateau
      print (measured_plateaus, mean_plateau, nominal_plateau )
      axs[2].hist ( measured_plateaus  , bins=100, range=( min([ mean_plateau,nominal_plateau] )*0.7 , max( [mean_plateau,nominal_plateau] ) *1.4    ), label='measured')
      axs[2].hist ( real_plateaus , bins=100 , label='real')
      aa = patches.Rectangle( (nominal_plateau_low, axs[0].get_ylim()[0]  ) , nominal_plateau_wid , axs[0].get_ylim()[1] , edgecolor='red' , facecolor='green' , alpha=0.2)
      axs[2].set_title("tidal  plateau [cc], measured: %2.1f [cc]"%(nominal_plateau))
      axs[2].legend(loc='upper left')
      axs[2].add_patch(aa)

      nominal_volume     =   np.nanmean( real_tidal_volumes)   # float ( meta[objname]["Tidal Volume"] ) / 10.
      print (nominal_volume)
      nominal_volume_low = nominal_volume - 4 - 0.15 * nominal_volume
      nominal_volume_wid = 8 + 0.3 * nominal_volume
      axs[3].hist ( measured_volumes  , bins=100, range=( min([ mean_volume,nominal_volume] )*0.7 , max( [mean_volume,nominal_volume] ) *1.4    ), label='measured')
      axs[3].hist ( real_tidal_volumes , bins=100 , label='real')
      aa = patches.Rectangle( (nominal_volume_low, axs[0].get_ylim()[0]  ) , nominal_volume_wid , axs[0].get_ylim()[1] , edgecolor='red' , facecolor='green' , alpha=0.2)
      axs[3].set_title("tidal  volume [cc], measured: %2.1f [cc],\nnominal %i [cc]"%(nominal_volume,int ( meta[objname]['Tidal Volume'])/10))
      axs[3].legend(loc='upper left')
      axs[3].add_patch(aa)

      figpath = "%s/%s_summary_%s.pdf" % (output_directory, meta[objname]['Campaign'], objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      figs.savefig(figpath)

      meta[objname]["mean_peep"]         =  mean_peep
      meta[objname]["rms_peep"]          =  rms_peep
      meta[objname]["mean_plateau"]      =  mean_plateau
      meta[objname]["rms_plateau"]       =  rms_plateau
      meta[objname]["mean_peak"]         =  mean_peak
      meta[objname]["rms_peak"]          =  rms_peak
      meta[objname]["mean_volume"]       =  mean_volume
      meta[objname]["rms_volume"]        =  rms_volume
      meta[objname]["simulator_volume"]  =  nominal_volume
      meta[objname]["simulator_plateau"] =  nominal_plateau

      filepath = "%s/summary_%s_%s.json" % (output_directory, meta[objname]['Campaign'],objname.replace('.txt', '')) # TODO: make sure it is correct, or will overwrite!
      json.dump( meta[objname], open(filepath , 'w' ) )


if __name__ == '__main__':
  import argparse
  import matplotlib
  import style

  parser = argparse.ArgumentParser(description='repack data taken in continuous mode')
  parser.add_argument("input", help="name of the MVM input file (.txt)")
  parser.add_argument("-d", "--output-directory", type=str, help="name of the output directory for plots", default="plots_iso")
  parser.add_argument("-i", "--ignore_sim", action='store_true',  help="ignore_sim")
  parser.add_argument("-skip", "--skip_files", type=str,  help="skip files", nargs='+', default="")
  parser.add_argument("-p", "--plot", action='store_true', help="show plots")
  parser.add_argument("-s", "--save", action='store_true', help="save HDF")
  parser.add_argument("-f", "--filename", type=str, help="single file to be processed", default='.')
  parser.add_argument("-c", "--campaign", type=str, help="single campaign to be processed", default="")
  parser.add_argument("-o", "--offset", type=float, help="offset between vent/sim", default='.0')
  parser.add_argument("--db-google-id", type=str, help="name of the Google spreadsheet ID for metadata", default="1aQjGTREc9e7ScwrTQEqHD2gmRy9LhDiVatWznZJdlqM")
  parser.add_argument("--db-range-name", type=str, help="name of the Google spreadsheet range for metadata", default="20200408 ISO!A2:AZ")
  parser.add_argument("--mvm-sep", type=str, help="separator between datetime and the rest in the MVM filename", default=" -> ")
  args = parser.parse_args()

  columns_rwa = ['dt',
    'airway_pressure',
    'muscle_pressure',
    'tracheal_pressure',
    'chamber1_vol',
    'chamber2_vol',
    'total_vol',
    'chamber1_pressure',
    'chamber2_pressure',
    'breath_fileno',
    'aux1',
    'aux2',
    'oxygen'
  ]
  columns_dta = [#'dt',
    'breath_no',
    'compressed_vol',
    'airway_pressure',
    'muscle_pressure',
    'total_vol',
    'total_flow',
    'chamber1_pressure',
    'chamber2_pressure',
    'chamber1_vol',
    'chamber2_vol',
    'chamber1_flow',
    'chamber2_flow',
    'tracheal_pressure',
    'ventilator_vol',
    'ventilator_flow',
    'ventilator_pressure',
  ]

  print ("INPUT : : ",  args.input )

  # read metadata spreadsheet
  df_spreadsheet = read_online_spreadsheet(spreadsheet_id=args.db_google_id, range_name=args.db_range_name)

  #if option -n, select only one test
  if len ( args.filename )  > 2 :
    unreduced_filename = args.filename.split("/")[-1]
    reduced_filename = '.'.join(unreduced_filename.split('.')[:-1])

    print ( "Selecting only: " ,  reduced_filename  )
    df_spreadsheet = df_spreadsheet[ ( df_spreadsheet["MVM_filename"] == unreduced_filename )  ]

  filenames = df_spreadsheet['MVM_filename'].unique()

  ntests = 0

  for filename in filenames:
    # continue if there is no filename
    if not filename: continue

    # read the metadata and create a dictionary with relevant info
    meta  = read_meta_from_spreadsheet (df_spreadsheet, filename)
    ntests += len(meta)

    objname = f'{filename}_0'   #at least first element is always there

    # compute the file location: local folder to the data repository + compaign folder + filename
    fname = f'{args.input}/{meta[objname]["Campaign"]}/{meta[objname]["MVM_filename"]}'
    if not fname.endswith(".txt"):
      fname = f'{fname}.txt'

    print(f'\nFile name {fname}')
    if fname.split('/')[-1] in args.skip_files:
      print('    ... skipped')
      continue

    if args.campaign:
      if args.campaign not in fname:
        print(f'    ... not in selected campaign {args.campaign}')
        continue

    # determine RWA and DTA data locations
    fullpath_rwa = f'{args.input}/{meta[objname]["Campaign"]}/{meta[objname]["SimulatorFileName"]}'

    if fullpath_rwa.endswith('.dta'):
      fullpath_rwa =  fullpath_rwa[:-4]      #remove extension if dta
    if not fullpath_rwa.endswith('.rwa'):
      fullpath_rwa =  f'{fullpath_rwa}.rwa'  #if .rwa extension not present, add it

    fullpath_dta = fullpath_rwa.replace('rwa', 'dta')
    print(f'will retrieve RWA and DTA simulator data from {fullpath_rwa} and {fullpath_dta}')

    # run
    process_run(meta, objname=objname, input_mvm=fname, fullpath_rwa=fullpath_rwa, fullpath_dta=fullpath_dta, columns_rwa=columns_rwa, columns_dta=columns_dta, save=args.save, manual_offset=args.offset,  ignore_sim=args.ignore_sim, mvm_sep=args.mvm_sep, output_directory=args.output_directory)

  if args.plot:
    if ( len (filenames) < 2 ) :
      plt.show()
    else :
      answer = input("plot all the files? (return: yes, Ctrl-D: no)")
      plt.show()
