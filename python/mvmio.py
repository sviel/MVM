import numpy as np
import pandas as pd

def get_raw_df(fname, columns, columns_to_deriv, timecol='dt'):
  df = pd.read_csv(fname, skiprows=4, names=columns, sep='\t', engine='python')
  for column in columns_to_deriv:
    df[f'deriv_{column}'] = df[column].diff() / df[timecol].diff() * 60. # TODO
  return df

def get_simulator_df(fullpath_rwa, fullpath_dta, columns_rwa, columns_dta):
  df_rwa = get_raw_df(fullpath_rwa, columns=columns_rwa, columns_to_deriv=['total_vol'])
  df_dta = get_raw_df(fullpath_dta, columns=columns_dta, columns_to_deriv=[])
  df0 = df_dta.join(df_rwa['dt'])
  df_rwa['oxygen'] = df_rwa['oxygen']  /  df_rwa['airway_pressure']
  df  = df0.join(df_rwa['oxygen'] )
  return df

def get_mvm_df(fname, sep=' -> '): 
  #data from the ventilator
  data = []

  is_unix = False

  with open(fname) as f:
    lines = f.readlines()
    for iline, line in enumerate(lines):
      if iline == 0: continue # skip first line

      line = line.strip()
      line = line.strip('\n')
      if not line: continue
      l = line.split(sep)
      try:
        par = sep.join(l[1:]).split(',')
      except :
        print (line)
        continue
      # remove unformatted lines
      try:
        for x in par: float(x)
      except ValueError:
        continue

      t = l[0]
      if ':' not in l[0]:
        t = float(l[0]) # in this way, t is either a string (if HHMMSS) or a float
        is_unix = True
      data.append({'date': t, 'flux':float(par[0]),'pressure':float(par[1]), 'airway_pressure':float(par[2]), 'in':float(par[3]),'flow2nd_der':float(par[4]), 'out':float(par[5])})

  df = pd.DataFrame(data)
  if not is_unix: # text timestamp
    df['dt'] = ( pd.to_datetime(df['date']) - pd.to_datetime(df['date'][0]) )/np.timedelta64(1,'s')
  else: # unix timestamp in seconds
    df['dt'] = ( pd.to_datetime(df['date'], unit='s') - pd.to_datetime(df['date'][0], unit='s') )/np.timedelta64(1,'s')
  #dtmax = df.iloc[-1,:]['dt']
  #timestamp = np.linspace( df.iloc[0,:]['dt'] ,  df.iloc[-1,:]['dt']*(dtmax-0.08)/dtmax , len(df) )   #use this line if you want to stretch the x axis of MVM data

  return df
