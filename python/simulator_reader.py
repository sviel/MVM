''' A simple reader and visualiser for RWA/DTA files '''
''' Valerio Ippolito - INFN Sezione di Roma '''

import pandas as pd
import matplotlib.pyplot as plt

def get_raw_df(fname, columns, columns_to_deriv, timecol='dt'):
  df = pd.read_csv(fname, skiprows=4, names=columns, sep='\t', engine='python')
  for column in columns_to_deriv:
    df[f'deriv_{column}'] = df[column].diff() / df[timecol].diff() * 60. # TODO
  return df

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser(description='repack data taken in continuous mode')
  parser.add_argument("input", help="name of the input file (.txt)", nargs='+')
# parser.add_argument("-o", "--output", type=str, help="name of the output folder", required=True)
# parser.add_argument("-d", "--duration", type=int, help="minimal duration of period aka subrun, in s", required=True)
  parser.add_argument("-p", "--plot", action='store_true', help="show plots")
  parser.add_argument("-s", "--save", action='store_true', help="save HDF")
  parser.add_argument("-f", "--folder", type=str, help="path to input data", default='.')
  parser.add_argument("-r", "--rwa-only", action='store_true', help="use RWA data only (when DTA is not available)")
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
  for fname in args.input:
    objname = f'{fname.split("/")[-1].replace(".txt", "")}'
    fullpath_rwa = f'{args.folder}/{fname}'
    fullpath_dta = fullpath_rwa.replace('rwa', 'dta')
    print(objname, fullpath_rwa, fullpath_dta)
    df_rwa = get_raw_df(fullpath_rwa, columns=columns_rwa, columns_to_deriv=['total_vol'])
    if not args.rwa_only:
      df_dta = get_raw_df(fullpath_dta, columns=columns_dta, columns_to_deriv=[])
      print(df_rwa.head())
      print(df_dta.head())
      df = df_dta.join(df_rwa['dt'])
    else:
      df = df_rwa
    print(df.head())
    
    if args.save:
      df.to_hdf(f'{objname}.h5', key='simulator')

    if args.plot:
      ax = df.plot(x='dt', y='muscle_pressure', label='muscle_pressure [cmH2O]')
      df.plot(ax=ax, x='dt', y='airway_pressure', label='airway_pressure [cmH2O]')
#     df.plot(ax=ax, x='dt', y='tracheal_pressure', label='tracheal_pressure')
     #df.plot(ax=ax, x='dt', y='total_vol', label='total_vol')
      plt.plot(df['dt'], df['total_vol']/10., label='total_vol [cl]')
      if not args.rwa_only:
        df.plot(ax=ax, x='dt', y='total_flow', label='total_flow [l/min]')
      else:
        df.plot(ax=ax, x='dt', y='deriv_total_vol', label='deriv_total_vol [l/min]')
#     df.plot(ax=ax, x='dt', y='aux1', label='aux1')
#     df.plot(ax=ax, x='dt', y='aux2', label='aux2')
      plt.gcf().suptitle(objname)
      ax.legend(loc='upper left')
      plt.show()
