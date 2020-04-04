import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import colors
import csv
from scipy.interpolate import interp1d
import matplotlib.patches as patches

def read_meta_csv(fname, simformat='.rwa'):
  ''' get a dict containing run metadata '''
  meta = {}
  with open(fname) as f:
    reader = csv.DictReader(f)

    for line in reader:
      #meta[f'run_{line["HardwareRun"]}'] = {
      meta[line['VentilatorFile']] = {
        'Rate respiratio': line["Frequency"],
        'Ratio respiratio': '1:2',
        'Peep': line["PEEP"],
        'Cycles': line["Cycles"],
        'Pinspiratia': line["Pinspiration"],
        'SimulatorFileName': f'{line["SimulatorFile"]}{simformat}',
    }
  return meta

def read_mhra_csv(fname):
  dfmhra = pd.read_csv(fname, sep=',')
  dfmhra['plot'] = 0
  return dfmhra

def get_raw_df(fname, columns, columns_to_deriv, timecol='dt'):
  df = pd.read_csv(fname, skiprows=4, names=columns, sep='\t', engine='python')
  for column in columns_to_deriv:
    df[f'deriv_{column}'] = df[column].diff() / df[timecol].diff() * 60. # TODO
  return df

def process_run(meta, objname,input_mvm, fullpath_rwa, fullpath_dta, columns_rwa, columns_dta, offset=0., save=False, ignore_sim=False, mhracsv=None):
  df_rwa = get_raw_df(fullpath_rwa, columns=columns_rwa, columns_to_deriv=['total_vol'])
  df_dta = get_raw_df(fullpath_dta, columns=columns_dta, columns_to_deriv=[])
  df0 = df_dta.join(df_rwa['dt'])
  df_rwa['oxygen'] = df_rwa['oxygen']  /  df_rwa['airway_pressure']
  df  = df0.join(df_rwa['oxygen'] )

  print (df.head())
  if ignore_sim == True :  print ("AAAA")
  #data from the ventilator (run number hardcoded for now)
  data = []

  with open(input_mvm) as f:
    lines = f.readlines()
    for line in lines:
      line = line.strip()
      line = line.strip('\n')
      if not line: continue
      l = line.split(' -> ')
      try:
        par = l[1].split(',')
      except :
        print (line)
        continue
      # remove unformatted lines
      try:
        for x in par: float(x)
      except ValueError:
        continue
      data.append({'date':l[0], 'flux':float(par[0]),'pressure':float(par[1]), 'airway_pressure':float(par[2]), 'in':float(par[3]),'flow2nd_der':float(par[4]), 'out':float(par[5])})

  dfhd = pd.DataFrame(data)
  dfhd['dt'] = ( pd.to_datetime(dfhd['date']) - pd.to_datetime(dfhd['date'][0]) )/np.timedelta64(1,'s')
  dtmax = dfhd.iloc[-1,:]['dt']
  #timestamp = np.linspace( dfhd.iloc[0,:]['dt'] ,  dfhd.iloc[-1,:]['dt']*(dtmax-0.08)/dtmax , len(dfhd) )
  timestamp = np.linspace( dfhd.iloc[0,:]['dt'] ,  dfhd.iloc[-1,:]['dt'] , len(dfhd) )
  deltaT = timestamp[2] - timestamp[1]
  #correct for miscalibration of pressure sensor
  dfhd['pressure'] =dfhd['pressure']

  #correct for volume units (ml -> cl)
  df['total_vol'] = df['total_vol'] * 0.1
  #dfhd = dfhd[(dfhd['dt']>df['dt'].iloc[0]) & (dfhd['dt']<df['dt'].iloc[-1]) ]

  #set MVM flux to 0 when out valve is open
  dfhd['flux'] = np.where ( dfhd['out'] >50 , 0 , dfhd['flux'] )
  dfhd['flux'] = np.where ( dfhd['flux'] < 0 , 0 , dfhd['flux'] )

  ##################################
  #rough shift for plotting purposes only
  ##################################
  imax1 = dfhd[ ( dfhd['dt']<14 ) ] ['flux'].idxmax()
  tmax1 = dfhd['dt'].iloc[imax1]
  imax2 = df[  (df['dt']<7) & (df['dt']>0.2)   ] ['total_flow'].idxmax()
  tmax2 = df['dt'].iloc[imax2]
  shift  = tmax2 - tmax1
  print (tmax1, tmax2, shift)
  dfhd['dt'] = timestamp + shift  + offset # add arbitrary shift to better match data

  ##################################
  #true start time of breaths based on muscle_pressure (first negative time)
  ##################################
  #negative_times = df[ ( df['muscle_pressure']<0 ) ]['dt']    #array of negative times
  #start_times = [ float (negative_times.iloc[i])  for i in range(0,len(negative_times)) if negative_times.iloc[i]-negative_times.iloc[i-1] > 1.  or i == 0  ]   #select only the first negative per bunch

  ##################################
  #true start time of breaths based on muscle_pressure (first negative time)
  ##################################
  negative_times = dfhd[ ( dfhd['out']>50 ) ]['dt']    #array of negative times
  start_times = [ float (negative_times.iloc[i])  for i in range(0,len(negative_times)-1) if negative_times.iloc[i+1]-negative_times.iloc[i] > 0.1  or i == 0  ]   #select only the first negative per bunch
  print(start_times)

  ##################################
  #reaction time of system based on positive flow
  ##################################
  positive_flow_times = df[ ( df['total_flow']>5 ) ]['dt']   # constant threshold to avoid artifact
  #flow_start_time = [ float (positive_flow_times.iloc[i])  for i in range(1,len(positive_flow_times)) if positive_flow_times.iloc[i]-positive_flow_times.iloc[i-1] > 1. ]   #select only the first negative per bunch
  flow_start_time = []
  c = []
  for st in start_times :
    isFound = False
    for pf in positive_flow_times :
      if pf > st and isFound == False:
        isFound = True
        flow_start_time.append(pf)
    if isFound == False :
      flow_start_time.append(0)


  print (len ( flow_start_time) , len (flow_start_time) )
  reaction_times = np.array ([100 * ( ft - st ) for ft,st in zip (flow_start_time,start_times) ])
  reaction_times = np.where ( abs(reaction_times)>1000, 0,reaction_times  )

  """
  #volume calculation
  dfhd['volume']  = 0
  #np.where ( dfhd['out'] <50 , 0 , dfhd['flux'] )
  #for i in range ( 0, len ( start_times)-1 ) :
  for i in range ( len ( start_times ) - 1 ) :
    dftmp = dfhd[ ( dfhd['dt']>start_times[i] ) & ( dfhd['dt'] < start_times[i+1] )  ]
    print (len(dftmp) )
    for j in range ( 1, len (dftmp)) :
      deltat =  dftmp['dt'].iloc[j] - dftmp['dt'].iloc[j-1]
      mydt   =  dftmp['dt'].iloc[j]
      dfhd.loc[ dfhd.dt == mydt, 'volume'] = dftmp['flux'].iloc[j-1] *  deltat * 100 #

  """


  #print (reaction_time)
  ##################################
  #add time stamp of breath start to each sample
  ##################################
  df['start'] = 0
  df['tbreath'] = 0
  df['reaction_time'] = 0
  for s,f in zip (start_times,reaction_times) :
    df.loc[df.dt>s,'start']   = s
    df.loc[df.dt>s,'tbreath'] = df.dt - s
    df.loc[df.dt>s,'reaction_time'] = f

  dfhd['start'] = 0
  for s in start_times :
    dfhd.loc[dfhd.dt>s,'start']   = s

  """
  ##################################
  #select first period
  ##################################
  dfcycle_all =  df[(df['start']< 52) & (df['start']> 4.5 )]

  dfcycle = []
  cola =  ( cm.Oranges(np.linspace(0,1, 12)))
  colw =  ( cm.YlGn(np.linspace(0,1, 12)))
  for i in range (0,12) :
    dfcycle.append ( df[(df['start']<8.5+4.*i) & (df['start']> 4.5+ 4.*i )] )
  """

  ##################################
  # chunks
  ##################################

  cycle = df.groupby(['start']).mean()

  df['mean_pressure'] = 0

  for i,r in cycle.iterrows():
    df.loc[df.start == i, 'mean_pressure'] = r.total_flow

  df['norm_pressure']  = df['airway_pressure'] - df['mean_pressure']

  df['max_pressure'] = 0
  cycle = df.groupby(['start']).max()
  for i,r in cycle.iterrows():
    df.loc[df.start == i, 'max_pressure'] = r.total_flow

  ##################################
  # tidal volume
  ##################################

  dfhd['out_times'] = dfhd['out'].diff()
  dfhd['out_times'] = np.where(dfhd['out_times']>0, 1,0)
  volumes_times = dfhd[dfhd['out_times'] ==1]['dt'].unique()

  if dfhd['out'].iloc[0] > 0 :
    volumes_times2 = [0]
    for x in volumes_times: volumes_times2.append(x)
    volumes_times = volumes_times2
    print ("AAAADSASASASA")
  dfhd['out_times'] = 0
  for s,v in zip(start_times, volumes_times):
    print (s,v)
    dfhd.loc[(dfhd.dt>s) & (dfhd.dt<v),'out_times'] = 1

  dfhd['tidal_volume'] = 0
  for c in dfhd['start'].unique():
    cum = dfhd[dfhd['start']==c]['flux'].cumsum()
    dfhd.loc[dfhd.start == c, 'tidal_volume'] = cum

  dfhd['tidal_volume'] *= deltaT/60.*100
  dfhd['tidal_volume'] *= dfhd['out_times']

  ##################################
  # find runs
  ##################################

  df['mean_max_pressure'] = df['max_pressure'].rolling(4).mean()

  df['run'] = abs(df['mean_max_pressure'].diff())
  df['run'] = np.where(df['run'] > 2, 100, 0)

  starts  = df[df.run > 0]['dt'].values

  # add run number
  df['run'] = -1
  cc = 1
  for i in range(0,len(starts)-1):
    if starts[i+1] - starts[i] > 25:
      df.loc[(df['dt'] >= starts[i]) & (df['dt'] < starts[i+1]), 'run'] = cc
      cc += 1

  df['run'] = df['run']*10

  if save:
    df.to_hdf(f'{objname}.sim.h5', key='simulator')
    dfhd.to_hdf(f'{objname}.mvm.h5', key='simulator')
    #df[(df['dt']>410)&(df['dt']<427)].to_csv(f'{objname}.sim.csv')
    #dfhd[(dfhd['dt']>410)&(dfhd['dt']<427)].to_csv(f'{objname}.mvm.csv')


  if args.plot:
    ttitles = ["R = 5; C = 50","R = 20; C = 20", "R = 10; C = 50"," "]
    resistances = [5,20,10]
    compliances = [50,20,50]

    colors = {  "muscle_pressure": "#009933"  , #green
      "sim_airway_pressure": "#cc3300" ,# red
      "total_flow":"#ffb84d" , #orange
      "total_vol":"#009933" , #orange
      "reaction_time" : "#999999", #
      "pressure" : "black" , #  blue
      "vent_airway_pressure": "#003399" ,# blue
      "flux" : "#3399ff" #light blue
    }
    linw = 2

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
#   df.plot(ax=ax, x='dt', y='deriv_total_vol', label='deriv_total_vol [l/min]')
#   df.plot(ax=ax, x='dt', y='aux1', label='aux1')
#   df.plot(ax=ax, x='dt', y='aux2', label='aux2')
    #plt.gcf().suptitle(objname)
    #plt.plot(start_times, [0]*len(start_times), 'bo', label='real cycle start time')
    #df.plot(ax=ax, x='dt', y='reaction_time', label='reaction time [10 ms]', c=colors['reaction_time'])
    #plt.plot(   ,  100 *reaction_times,      label='reaction time ', marker='o', markersize=1, linewidth=0, c='red')
    #ax.
    for i,t in enumerate(start_times) :
      ax.text(t, 0.5, "%i"%i, verticalalignment='bottom', horizontalalignment='center', color='red', fontsize=35)

    ax.set_xlabel("Time [sec]")
    ax.legend(loc='upper center', ncol=2)

    """
    #a clean canavs with simulator only data
    #fig5,ax5 = plt.subplots()
    ax5 = df.plot(x='dt', y='airway_pressure', label='airway_pressure [cmH2O]', c=colors['sim_airway_pressure'], linewidth = linw)
    #df.plot(ax=ax5, x='dt', y='muscle_pressure', label='muscle_pressure [cmH2O]', c=colors['muscle_pressure'], linewidth = linw)
    df.plot(ax=ax5, x='dt', y='run', label='run index')
    df.plot(ax=ax5, x='dt', y='total_vol', label='total_vol [l]', c=colors['total_vol'], linewidth = linw)
    df.plot(ax=ax5, x='dt', y='total_flow', label='total_flow [l/min]', c=colors['total_flow'], linewidth = linw)
    #df.plot(ax=ax5, x='dt', y='oxygen', label='oxygen [%]', c='red', linewidth = linw)
    plt.gcf().suptitle("simulator only data")
    plt.plot(start_times, [0]*len(start_times), 'bo', label='real cycle start time')
    ax5.legend(loc='lower left')

    #a clean canavs with ventilator only data
    #fig5b,ax5b = plt.subplots()
    ax5b = dfhd.plot( x='dt', y='pressure', label='ventilator pressure [cmH2O]', c=colors['pressure'])
    dfhd.plot(ax=ax5b, x='dt', y='airway_pressure', label='ventilator airway pressure [cmH2O]', c=colors['vent_airway_pressure'])
    dfhd.plot(ax=ax5b, x='dt', y='flux', label='ventilator flux [l/min]', c=colors['flux'])
    #dfhd.plot(ax=ax5b, x='dt', y='in', label='in', c='green')
    #dfhd.plot(ax=ax5b, x='dt', y='flow2nd_der', label='flow2nd_der?', c='pink')
    #dfhd.plot(ax=ax5b, x='dt', y='out', label='out', c='red')
    plt.gcf().suptitle("ventilator only data")
    ax5b.legend(loc='lower left')


    maxrunidx = len( df['run'].unique() ) - 1
    fig8,ax8 = plt.subplots(2,2)
    pad8 = ax8.flatten ()
    for i in range (0, maxrunidx) :
      #make a subset dataframe for simulator
      dftmp = df[(df.run==(i+1)*10)]                    #select corresponding period
      start_times_in_run = dftmp['start'].unique()      #array of start times in period
      #select three central cycles in period
      dftmp = dftmp[(dftmp.start==start_times_in_run[6])|(dftmp.start==start_times_in_run[7])|(dftmp.start==start_times_in_run[8]) ]
      dftmp.plot(ax=pad8[i], x='dt', y='airway_pressure', label='airway_pressure [cmH2O]', c=colors['sim_airway_pressure'], linewidth = linw)

      first_time_bin  = dftmp['dt'].iloc[0]
      last_time_bin   = dftmp['dt'].iloc[len(dftmp)-1]
      #make a subset dataframe for ventilator
      dfvent = dfhd[ (dfhd['dt']>first_time_bin) & (dfhd['dt']<last_time_bin) ]
      dfvent.plot(ax=pad8[i],  x='dt', y='pressure', label='ventilator pressure [cmH2O]', c=colors['pressure'], linewidth = linw)
      dfvent.plot(ax=pad8[i],  x='dt', y='airway_pressure', label='ventilator airway pressure [cmH2O]', c=colors['vent_airway_pressure'])
      dfvent.plot(ax=pad8[i],  x='dt', y='flux', label='ventilator flux [l/min]', c=colors['flux'] , linewidth = linw)
      dftmp.plot(ax=pad8[i], x='dt', y='total_flow', label='total_flow [l/min]', c=colors['total_flow'], linewidth = linw)
      pad8[i].set_title("%s; PEEP%s; P_i%s; R%s"%(ttitles[i], meta[objname]["Peep"], meta[objname]["Pinspiratia"], meta[objname]["Rate respiratio"]))
      #pad8[i].set_title('test')
    """

    three_start_times = [20,50,120]

    for i in range (0, 3) :
      skip_plot = False ;

      if mhracsv is not None :
        PE = meta[objname]["Peep"]
        PI = meta[objname]["Pinspiratia"]
        RR = meta[objname]["Rate respiratio"]
        print ("Looking for %s cycle, R=%s, C=%s, RR=%s, PEEP=%s, PINSP=%s"%(0,resistances[i],compliances[i],RR,PE,PI) )

        condition = "(( dfmhra['R']==resistances[i])  & (dfmhra['RR']==int(RR)) & (dfmhra['C']==compliances[i]) & (dfmhra['PEEP']==int ( PE) ) & (dfmhra['PINSP']==int ( PI))"

        my_selected_cycle_series =  mhracsv[( mhracsv['R']==resistances[i])  & (mhracsv['RR']==int(RR)) & (mhracsv['C']==compliances[i]) & (mhracsv['PEEP']==int ( PE) ) & (mhracsv['PINSP']==int ( PI)) ]['CYCLE_START']
        print (len ( my_selected_cycle_series))

        if ( len( my_selected_cycle_series)==1 ) :
          my_selected_cycle = my_selected_cycle_series.iloc[0]
          print (my_selected_cycle)
          mhracsv.loc[ ( mhracsv['R']==resistances[i])  & (mhracsv['RR']==int(RR)) & (mhracsv['C']==compliances[i]) & (mhracsv['PEEP']==int ( PE) ) & (mhracsv['PINSP']==int ( PI))      ,    'plot'] = 1
          three_start_times = [my_selected_cycle, my_selected_cycle+1 ,my_selected_cycle+2]
        else : skip_plot = True

      if skip_plot == True : continue

      fig11,ax11 = plt.subplots()

      #make a subset dataframe for simulator
      dftmp = df[ df['start']> three_start_times[i] ]
      start_times_in_run = dftmp['start'].unique()      #array of start times in period
      dftmp = dftmp[(dftmp.start>start_times_in_run[0])&(dftmp.start<start_times_in_run[4])]
      #select three central cycles in period
      first_time_bin  = dftmp['dt'].iloc[0]
      last_time_bin   = dftmp['dt'].iloc[len(dftmp)-1]
      #make a subset dataframe for ventilator
      dfvent = dfhd[ (dfhd['dt']>first_time_bin) & (dfhd['dt']<last_time_bin) ]

      dfvent.plot(ax=ax11,  x='dt', y='tidal_volume',    label='MVM tidal volume       [cl]', c='pink')
      dfvent.plot(ax=ax11,  x='dt', y='airway_pressure', label='MVM airway pressure [cmH2O]', c=colors['vent_airway_pressure'])
      dfvent.plot(ax=ax11,  x='dt', y='flux',            label='MVM flux            [l/min]', c=colors['flux'])
      dftmp.plot(ax=ax11, x='dt', y='airway_pressure',   label='SIM airway pressure [cmH2O]', c=colors['sim_airway_pressure'])
      dftmp.plot(ax=ax11, x='dt', y='total_flow',        label='SIM flux            [l/min]', c=colors['total_flow'])
      #dftmp.plot(ax=ax11, x='dt', y='total_vol',         label='SIM volume             [dl]', c=colors['total_vol'] , linestyle='--')
      #dfvent.plot(ax=ax11,  x='dt', y='out',             label='MVM out                    ', c='green')
      #dfvent.plot(ax=ax11, x='dt', y='volume',           label='vol                [ml]', c='green')
      #ax11.set_title("%s; PEEP%s; P_i%s; R%s"%(ttitles[i], meta[objname]["Peep"], meta[objname]["Pinspiratia"], meta[objname]["Rate respiratio"]))
      ymin, ymax = ax11.get_ylim()
      ax11.set_ylim(ymin*1.4, ymax*1.5)
      ax11.legend(loc='upper center', ncol=2)
      title1="R = %i [cmH2O/l/s]       C = %i [ml/cmH20]       PEEP = %s [cmH20]"%(resistances[i],
        compliances[i],
        meta[objname]["Peep"]
      )
      title2="Inspiration Pressure = %s [cmH20]      Frequency = %s [breath/min]"%(
        meta[objname]["Pinspiratia"],
        meta[objname]["Rate respiratio"]
      )
      ax11.set_xlabel("Time [s]")
      ax11.text(0.5, 0.08, title2, verticalalignment='bottom', horizontalalignment='center', transform=ax.transAxes, color='#7697c4', fontsize=20)
      ax11.text(0.5, 0.026, title1, verticalalignment='bottom', horizontalalignment='center', transform=ax.transAxes, color='#7697c4', fontsize=20)
      nom_pressure = float(meta[objname]["Pinspiratia"])
      xmin, xmax = ax11.get_xlim()
      rect = patches.Rectangle((xmin,nom_pressure-2),xmax-xmin,4,edgecolor='None',facecolor='b', alpha=0.2)
      print ("nominal insp press", nom_pressure, xmin, xmax)


      fig11.savefig("plots/%s_%i.png"%(objname,i))

    """
    fig6,ax6 = plt.subplots(3,3)
    pad6 = ax6.flatten ()
    fig7,ax7 = plt.subplots(3,3)
    pad7 = ax7.flatten ()

    df['dt0']     = df['tbreath']
    #df['abs_dt0'] = 0
    #df['abs_dt0'] = [ pd.Timedelta( x , unit ='s'  ) for x in df['tbreath'] ]
    df['abs_dt0'] = pd.to_timedelta(df['dt0'],'s')

#    for i,r in df.iterrows():
#      df.iloc[i, 'abs_dt0'] = pd.Timedelta( r.dt0, unit ='s'  )
    #resamp = df.resample(rule='0.01S' , on='abs_dt0')

    for i in range (0, maxrunidx) :

      #bins = np.linespace(0,4.,1000)
      #mdftmp['tbreath_digitised'] = np.digitize(mdftmp['tbreath'], bins)
      dftmp = df[(df.run==i+1)].sort_values('tbreath')
      # remove first and last cycle from each run
      c = dftmp['start'].unique()
      dftmp.loc[dftmp['start']==c[0],  'run'] = -1
      dftmp.loc[dftmp['start']==c[-1], 'run'] = -1

      #mdftmp = dftmp[dftmp.run==i+1].groupby('tbreath').mean()

      dftmp['idx'] = [i for i,el in enumerate(dftmp['tbreath']) ]
      mdftmp = dftmp[dftmp.run==i+1].groupby('tbreath').mean()

      dftmp[dftmp.run==i+1].plot(ax=pad6[i], x='tbreath', y='airway_pressure', label='airway_pressure [cmH2O]', marker='.', markersize=0.15, linewidth=0, c=colors['airway_pressure'])
      dftmp[dftmp.run==i+1].plot(ax=pad6[i], x='tbreath', y='total_flow', label='total_flow [l/min]', marker='.', markersize=0.15, linewidth=0, c=colors['total_flow'] )
      #mdftmp.plot(ax=pad6[i], x='dt0', y='airway_pressure', label='total_flow [l/min]', c='black', marker='.', markersize=0.1, linewidth=0 )
      """

    """
    fig2,ax2 = plt.subplots()
    for i in range (0,12) :
      dfcycle[i].plot(ax=ax2, x='tbreath', y='total_flow',      label='total_flow [l/min]', marker='.', markersize=0.3, linewidth=0, c=cola[i])
      dfcycle[i].plot(ax=ax2, x='tbreath', y='airway_pressure', label='airway_pressure [cmH2O]', marker='.', markersize=0.2, linewidth=0 ,  c=colw[i])
      dfcycle[i].plot(ax=ax2, x='tbreath', y='muscle_pressure', label='muscle_pressure  [cmH2O]', marker='.', markersize='0.1', linewidth=0, c='blue')
    ax2.legend(loc='upper right')

    fig3,ax3 = plt.subplots()
    dfcycle_all.plot(ax=ax3, x='tbreath', y='total_flow',      label='total_flow [l/min]', marker='.', markersize=0.4, linewidth=0, c='blue')
    dfcycle_all.plot(ax=ax3, x='tbreath', y='airway_pressure', label='airway_pressure [cmH2O]', marker='.', markersize=0.4, linewidth=0 ,  c='red')
    dfcycle_all.plot(ax=ax3, x='tbreath', y='muscle_pressure', label='muscle_pressure  [cmH2O]', marker='.', markersize='0.4', linewidth=0, c='green')
    ax3.legend(loc='upper right')
    """
    #fig4,ax4 = plt.subplots()
    #plt.plot(start_times,  reaction_times,      label='tdiff', marker='o', markersize=1, linewidth=0, c='red')
    #ax4.legend(loc='upper right')


    plt.show()


if __name__ == '__main__':
  import argparse
  import matplotlib
  import style
  #plt.subplots(dpi=72)




  matplotlib.rcParams['font.sans-serif'] = "Courier New"
  # Then, "ALWAYS use sans-serif fonts"
  matplotlib.rcParams['font.family'] = "sans-serif"


  parser = argparse.ArgumentParser(description='repack data taken in continuous mode')
  parser.add_argument("input", help="name of the MVM input file (.txt)", nargs='+')
# parser.add_argument("-o", "--output", type=str, help="name of the output folder", required=True)
# parser.add_argument("-d", "--duration", type=int, help="minimal duration of period aka subrun, in s", required=True)
  parser.add_argument("-i", "--ignore_sim", action='store_true',  help="ignore_sim")
  parser.add_argument("-p", "--plot", action='store_true', help="show plots")
  parser.add_argument("-s", "--save", action='store_true', help="save HDF")
  parser.add_argument("-m", "--mhra", action='store_true', help="save mhra")
  parser.add_argument("-f", "--folder", type=str, help="path to input simulator data", default='.')
  parser.add_argument("-o", "--offset", type=float, help="offset between vent/sim", default='.0')
  parser.add_argument("-l", "--logbook", type=str, help="logbook location", default='../Data/logbook.csv')
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

  # read logbook with association SIM-MVM
  meta = read_meta_csv(args.logbook)
  mhra = None
  if (args.mhra==True) : mhra = read_mhra_csv("../Data/logbook.MHRA.csv")

  for fname in args.input:
    # retrieve run name
    objname = f'{fname.split("/")[-1]}'
    print (fname, objname )
    print(f'processing objname={objname}')
    fullpath_rwa = f'{args.folder}/{meta[objname]["SimulatorFileName"]}'
    fullpath_dta = fullpath_rwa.replace('rwa', 'dta')
    print(f'will retrieve RWA and DTA simualtor data from {fullpath_rwa} and {fullpath_dta}')

    # run
    process_run(meta, objname=objname, input_mvm=fname, fullpath_rwa=fullpath_rwa, fullpath_dta=fullpath_dta, columns_rwa=columns_rwa, columns_dta=columns_dta, save=args.save, offset=args.offset,  ignore_sim=args.ignore_sim, mhracsv=mhra)

  if args.mhra==True :
    mhra.to_csv("mhra.out.csv")
