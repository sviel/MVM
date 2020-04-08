''' Read JSON output from combine.py and prepare summary plots '''

import pandas as pd
import json
import matplotlib.pyplot as plt
import lmfit

def process_files(files):
  dfs = []
  for i, fname in enumerate(files):
    dfs.append(pd.DataFrame(json.loads(open(fname).read()), index=[i]))
  df = pd.concat(dfs)

  df['Tidal Volume'] = df['Tidal Volume'].astype(int)
  df['simulator_volume_ml'] = df['simulator_volume'] * 10 # TODO patch conversion from ml to cl, remove when appropriate
  df['mean_volume_ml'] = df['mean_volume'] * 10
  df['rms_volume_ml'] = df['rms_volume'] * 10

 # Index(['Campaign', 'Compliance', 'I:E', 'MVM_filename', 'Peep', 'Pinspiratia',
 #      'Rate respiratio', 'Resistance', 'Run', 'SimulatorFileName',
 #      'Tidal Volume', 'cycle_index', 'mean_peak', 'mean_peep', 'mean_plateau',
 #      'mean_volume', 'rms_peak', 'rms_peep', 'rms_plateau', 'rms_volume',
 #      'simulator_volume', 'test_name'],

  TV = {
    '$V_{tidal}$ > 300 ml': df[df['Tidal Volume'] > 300],
    '300 ml > $V_{tidal}$ > 50 ml': df[(300 > df['Tidal Volume']) & (df['Tidal Volume'] > 50)],
    '$V_{tidal}$ < 50 ml': df[df['Tidal Volume'] < 50]
  }

  variables = [
    ('set BAP [$cmH_{2}O$]', 'measured BAP [$cmH_{2}O$]', 'Peep', 'mean_peep', 'rms_peep'),
    ('set $P_{insp}$ [$cmH_{2}O$]', 'measured $P_{plateau}$ [$cmH_{2}O$]', 'Pinspiratia', 'mean_plateau', 'rms_plateau'),
    ('set $P_{insp}$ [$cmH_{2}O$]', 'measured $P_{peak}$ [$cmH_{2}O$]', 'Pinspiratia', 'mean_peak', 'rms_peak'),
    ('set TV [ml]', 'measured TV [ml]', 'Tidal Volume', 'mean_volume_ml', 'rms_volume_ml'),
    ('TV from simulator [ml]', 'measured TV [ml]', 'simulator_volume_ml', 'mean_volume_ml', 'rms_volume_ml'),
  ]
  
  line = lmfit.models.LinearModel()
  for xname, yname, setval, mean, rms in variables:
    fig, ax = plt.subplots(1, 1)
    fig.canvas.set_window_title(f'x={setval}, y={mean}, yerr={rms}')

    for setname, data in TV.items():
#     params = line.guess(data[mean], x=data[setval])
#     res = line.fit(data[mean], params, x=data[setval], weights=1./data[rms])

      ax.errorbar(data[setval], data[mean], yerr=data[rms], fmt='o', label=setname)

    # linear fit
    df_to_fit = df
    if setval == 'simulator_volume_ml':
      df_to_fit = df_to_fit[df_to_fit[setval] > 50] # 201.12.1.104 from ISO
    params = line.guess(df_to_fit[mean], x=df_to_fit[setval])
    res = line.fit(df_to_fit[mean], params, x=df_to_fit[setval], weights=1./df_to_fit[rms])

    print(res.fit_report())
   #fitstring = f'${res.best_values["intercept"]} \pm {(res.best_values["slope"]-1)*100}$%'
    fitstring = f'$\pm$({res.best_values["intercept"]:.1f} +({(res.best_values["slope"]-1)*100:.0f}% of the value))'
    ax.plot(df_to_fit[setval], res.best_fit, '-', label=fitstring)
    ax.legend()
    ax.set_xlabel(f'{xname}')
    ax.set_ylabel(f'{yname}')
    #fig.show()
  plt.show()
  

if __name__ == '__main__':
  import argparse
  import style

  parser = argparse.ArgumentParser(description='prepare run summary tables from JSON')
  parser.add_argument("input", help="name of the input file (.txt)", nargs='+')
  parser.add_argument("-p", "--plot", action='store_true', help="show plots")
  args = parser.parse_args()

  process_files(args.input)
