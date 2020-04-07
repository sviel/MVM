import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import colors
from scipy.interpolate import interp1d
import matplotlib.patches as patches

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
  df['pressure'] = df['pressure']  + pressure_offset

  #set MVM flux to 0 when out valve is open AND when the flux is negative
  df['flux'] = np.where ( df['out'] > pv2_thr , 0 , df['flux'] )
  df['flux'] = np.where ( df['flux'] < 0 , 0 , df['flux'] )

def apply_rough_shift(sim, mvm, manual_offset):
  imax1 = mvm[ ( mvm['dt']<17 ) & (mvm['dt']>0.2) ] ['flux'].idxmax()
  tmax1 = mvm['dt'].iloc[imax1]
  imax2 = sim[  (sim['dt']<17) & (sim['dt']>0.2)   ] ['total_flow'].idxmax()
  tmax2 = sim['dt'].iloc[imax2]
  shift  = tmax2 - tmax1
  mvm['dt'] = mvm['timestamp'] + shift  + manual_offset # add arbitrary shift to better match data

def add_pv2_status(df):
  df['out_diff'] = df['out'].diff().fillna(0)
  df['out_status'] = np.where(df['out_diff'] == 0, 'steady', np.where(df['out_diff'] < 0, 'closing', 'opening'))

def get_start_times(df, thr=50, dist=0.1, quantity='out', timecol='dt'):
  ''' Get times where a given quantity starts being above threshold '''
  times_open = df[ ( df[quantity]>thr ) ][timecol]
  start_times = [ float(times_open.iloc[i]) for i in range(0, len(times_open)-1) if times_open.iloc[i+1]-times_open.iloc[i] > 0.1  or i == 0  ]
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
  for s in start_times :
    mvm.loc[mvm.dt>s,'start']   = s

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

def add_clinical_values (df) :
  """Add for reference the measurement of "TRUE" clinical values as measured using the simulator"""

  #true resistance
  df['delta_pout_pin']        = df['airway_pressure'] - df['chamber1_pressure']
  df['delta_vol']             = ( df['chamber1_vol'] + df['chamber2_vol'] ) .diff()
  df['airway_resistance']     = df['delta_pout_pin'] / df['delta_vol']
  df.loc[ (abs(df.airway_resistance)>250) | (df.airway_resistance<0) ,"airway_resistance"] = 0

  #true compliance
  df['deltapin']     =  df['chamber1_pressure'].diff()
  df['delta_volin']  =  ( df['chamber1_vol']  + df['chamber2_vol']) .diff()
  df['compliance']   =  df['delta_volin']/df['deltapin']
  df.loc[abs(df.compliance)>100,"compliance"] = 0


def measure_clinical_values(df, start_times):
  ''' Compute tidal volume and other clinical quantities for MVM data '''
  deltaT = get_deltat(df)

  # determine inspiration end times
  inspiration_end_times = df[df['out_status'] == 'opening']['dt'].unique()
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

  for s,v in zip(start_times, inspiration_end_times):
    this_inspiration = (df.dt>s) & (df.dt<v) # all samples
    first_sample = (df.dt == s) # start of inspiration
    last_sample  = (df.dt == v) # beginning of expiration

    df.loc[this_inspiration, 'is_inspiration'] = 1

#   print(df[this_inspiration]['timestamp'], df[this_inspiration]['dt'], df[this_inspiration]['pressure'].max())
    print(df[last_sample]['timestamp'], df[last_sample]['dt'], df[last_sample]['pressure'])
    df.loc[this_inspiration, 'cycle_tidal_volume'] = df[this_inspiration]['flux'].sum() * deltaT/60. * 100
    df.loc[this_inspiration, 'cycle_peak_pressure'] = df[this_inspiration]['pressure'].max()
    df.loc[this_inspiration, 'cycle_plateau_pressure'] = df[last_sample]['pressure']
    df.loc[this_inspiration, 'cycle_PEEP'] = df[first_sample]['pressure']

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

def process_run(meta, objname, input_mvm, fullpath_rwa, fullpath_dta, columns_rwa, columns_dta, manual_offset=0., save=False, ignore_sim=False, mhracsv=None, pressure_offset=0, mvm_sep=' -> '):
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


  ##################################
  # cycles
  ##################################

  # add PV2 status info
  add_pv2_status(dfhd)

  # compute cycle start
  # start_times = get_muscle_start_times(df) # based on muscle pressure
  start_times = get_start_times(dfhd) # based on PV2
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

  print('test:', dfhd['tidal_volume'].max())
  print('TV:', dfhd['cycle_tidal_volume'].unique(), dfhd['cycle_tidal_volume'].mean())
  print('Ppeak:', dfhd['cycle_peak_pressure'].unique(), dfhd['cycle_peak_pressure'].mean())
  print('Pplateau:', dfhd['cycle_plateau_pressure'].unique(), dfhd['cycle_plateau_pressure'].mean())
  print('PEEP:', dfhd['cycle_PEEP'].unique(), dfhd['cycle_PEEP'].mean())

  #plt.plot(dfhd['dt'], dfhd['out'], label='control out')
  plt.plot(dfhd['dt'], dfhd['pressure'], label='pressure')
  plt.plot(dfhd['dt'], dfhd['is_inspiration'], label='is_inspiration')
  #plt.plot(dfhd['dt'], dfhd['cycle_tidal_volume'], label='cycle_tidal_volume')
  #plt.plot(dfhd['dt'], dfhd['cycle_PEEP'], label='cycle_PEEP')
  plt.plot(dfhd['dt'], dfhd['tidal_volume'], label='tidal_volume')

  #plt.plot(df['dt'], df['tracheal_pressure'], label='true tracheal_pressure')
  #plt.plot(df['dt'], df['airway_pressure'], label='true airway_pressure')
  #plt.plot(df['dt'], df['chamber1_pressure'], label='true ch1_pressure')
  #plt.plot(df['dt'], df['chamber2_pressure'], label='true ch2_pressure')
  plt.plot(df['dt'], df['airway_resistance'], label='true airway_resistance')
  plt.plot(df['dt'], df['compliance'], label='true compliance')

  plt.legend()
  plt.show()
  print("Debugging ended - Press enter to continue to plots")
  input()
#  raise RuntimeError('Valerio')

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
    dfhd.plot(ax=ax, x='dt', y='out', label='out')
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
      skip_plot = False ;

      local_objname = "%s_%i"% ( objname[:-2] , i )

      PE = meta[local_objname]["Peep"]
      PI = meta[local_objname]["Pinspiratia"]
      RR = meta[local_objname]["Rate respiratio"]
      RT = meta[local_objname]["Resistance"]
      CM = meta[local_objname]["Compliance"]
      print ("Looking for R=%s, C=%s, RR=%s, PEEP=%s, PINSP=%s"%(RT,CM,RR,PE,PI) )

      my_selected_cycle = meta[local_objname]["cycle_index"]

      print ("\nFor test [ %s ]  I am selecting cycle %i, starting at %f \n"%(meta[local_objname]["test_name"], my_selected_cycle , start_times[ my_selected_cycle ]))

      fig11,ax11 = plt.subplots()

      #make a subset dataframe for simulator
      dftmp = df[ (df['start'] >= start_times[ my_selected_cycle ] ) & ( df['start'] < start_times[ my_selected_cycle + 3])  ]

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
      ax11.set_ylim(ymin*1.4, ymax*1.5)
      ax11.legend(loc='upper center', ncol=2)
      title1="R = %i [cmH2O/l/s]       C = %i [ml/cmH20]       PEEP = %s [cmH20]"%(RT,CM,PE )
      title2="Inspiration Pressure = %s [cmH20]      Frequency = %s [breath/min]"%(PI,RR)

      ax11.set_xlabel("Time [s]")

      xmin, xmax = ax11.get_xlim()
      ymin, ymax = ax11.get_ylim()
      ax11.text((xmax-xmin)/2.+xmin, 0.08*(ymax-ymin) + ymin,   title2, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      ax11.text((xmax-xmin)/2.+xmin, 0.026*(ymax-ymin) + ymin,  title1, verticalalignment='bottom', horizontalalignment='center', color='#7697c4')
      nom_pressure = float(meta[local_objname]["Pinspiratia"])
      rect = patches.Rectangle((xmin,nom_pressure-2),xmax-xmin,4,edgecolor='None',facecolor='green', alpha=0.2)
      ax11.add_patch(rect)

      fig11.savefig("plots_mhra/%s_C%i_R%i.pdf"%(local_objname, CM, RT))



if __name__ == '__main__':
  import argparse
  import matplotlib
  import style

  parser = argparse.ArgumentParser(description='repack data taken in continuous mode')
  parser.add_argument("input", help="name of the MVM input file (.txt)")
# parser.add_argument("-o", "--output", type=str, help="name of the output folder", required=True)
# parser.add_argument("-d", "--duration", type=int, help="minimal duration of period aka subrun, in s", required=True)
  parser.add_argument("-i", "--ignore_sim", action='store_true',  help="ignore_sim")
  parser.add_argument("-p", "--plot", action='store_true', help="show plots")
  parser.add_argument("-s", "--save", action='store_true', help="save HDF")
  parser.add_argument("-f", "--filename", type=str, help="single file to be processed", default='.')
  parser.add_argument("-o", "--offset", type=float, help="offset between vent/sim", default='.0')
  parser.add_argument("--db-google-id", type=str, help="name of the Google spreadsheet ID for metadata", default="1aQjGTREc9e7ScwrTQEqHD2gmRy9LhDiVatWznZJdlqM")
  parser.add_argument("--db-range-name", type=str, help="name of the Google spreadsheet range for metadata", default="20200407 ISO!A2:AZ")
  parser.add_argument("--mvm-sep", type=str, help="separator between datetime and the rest in the MVM filename", default=" -> ")
  args = parser.parse_args()

  columns_rwa = ['dt', 'airway_pressure', 'muscle_pressure', 'tracheal_pressure', 'chamber1_vol', 'chamber2_vol', 'total_vol', 'chamber1_pressure', 'chamber2_pressure', 'breath_fileno', 'aux1', 'aux2', 'oxygen']
  columns_dta = [#'dt',
    'breath_no',
    'compressed_vol',
    'airway_pressure',
    'muscle_pressure',
    'total_vol',
    'total_flow',
    'chamber1_pressure', 'chamber2_pressure',
    'chamber1_vol', 'chamber2_vol',
    'chamber1_flow', 'chamber2_flow',
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
    #continue if there is no filename
    if not filename: continue

    #read the metadata and create a dictionary with relevant info
    meta  = read_meta_from_spreadsheet (df_spreadsheet, filename)
    ntests += len(meta)

    objname = f'{filename}_0'   #at least first element is always there

    #compute the file location: local folder to the data repository + compaign folder + filename
    fname = f'{args.input}/{meta[objname]["Campaign"]}/{meta[objname]["MVM_filename"]}'
    if not fname.endswith(".txt"):
      fname = f'{fname}.txt'

    # determine RWA and DTA data locations
    fullpath_rwa = f'{args.input}/{meta[objname]["Campaign"]}/{meta[objname]["SimulatorFileName"]}'
    if fullpath_rwa.endswith('.dta'):
      fullpath_rwa =  fullpath_rwa[:-4]      #remove extension if dta
    if not fullpath_rwa.endswith('.rwa'):
      fullpath_rwa =  f'{fullpath_rwa}.rwa'  #if .rwa extension not present, add it

    fullpath_dta = fullpath_rwa.replace('rwa', 'dta')
    print(f'will retrieve RWA and DTA simulator data from {fullpath_rwa} and {fullpath_dta}')

    # run
    process_run(meta, objname=objname, input_mvm=fname, fullpath_rwa=fullpath_rwa, fullpath_dta=fullpath_dta, columns_rwa=columns_rwa, columns_dta=columns_dta, save=args.save, manual_offset=args.offset,  ignore_sim=args.ignore_sim, mvm_sep=args.mvm_sep)
    if args.plot:
      if ( len (filenames) < 2 ) :
        plt.show()
      else :
        answer = input("plot all the files? (return: yes, Ctrl-D: no)")
        plt.show()
