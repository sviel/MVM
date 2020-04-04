''' Put plots together in LaTeX documents '''

import re
from combine import read_meta_csv, read_mhra_csv

def associate_plots(meta, mhra, files, resistances, compliances, regexp):
  # loop over available files
  mhra['fname'] = ''

  for fname in files:
    m = regexp.match(fname)
    if m:
      plotset = int(m.group('plotset'))
      this_meta = meta[m.group('datafile')]
      RR = int(this_meta['Rate respiratio'])
      PE = int(this_meta['Peep'])
      PI = int(this_meta['Pinspiratia'])
      sel = (mhra['R'] == resistances[plotset]) & (mhra['RR'] == RR) & (mhra['C'] == compliances[plotset]) & (mhra['PEEP'] == PE) & (mhra['PINSP'] == PI)
      
      mhra.loc[sel, 'fname'] = fname
      mhra.loc[sel, 'plot'] = 1
    else:
      raise RuntimeError(f'Unknown MVM file {fname}')

  return mhra

def print_latex(mhra, title, rows=3, columns=2):
  toprint = []

  toprint.append(r'''
    \documentclass[a4paper]{article}
    \usepackage{subcaption}
    \usepackage{graphicx}
    \usepackage{todonotes}
    \usepackage{grffile} % for graphics extension
    \usepackage{siunitx} % for units
    \usepackage[margin=1.0cm]{geometry}
  ''')
  toprint.append(f'\\title{{{title}}}\n')
  toprint.append(r'''
    \begin{document}
    \maketitle
  ''')

  for i, row in mhra.iterrows():
    newfig = i % (rows * columns) == 0
    fname_form = row["fname"].split('/')[-1] # relative path to current dir

    if newfig:
      if i != 0:
        toprint.append(r'''  \caption{I:E=1:2, $O_2$=90-100}
          \end{figure}
        ''')
      toprint.append(r'''
      \begin{figure}[h!]
      \centering
      ''')

    if i % columns == 0 and not newfig:
      toprint.append(r'''
      ''')
    toprint.append(r''' \begin{subfigure}[b]{0.45\linewidth}'''
    )
    if row['plot'] == 1:
      toprint.append(r'''
        \includegraphics[width=\linewidth]'''
      )
      toprint.append(f'{{{fname_form}}}\n')
    else:
      toprint.append(r'''
        \missingfigure[figwidth=\linewidth]'''
      )
      toprint.append(f'{{Missing test}}\n')
    toprint.append(r'''      \caption{''')
    toprint.append(f'%%% $C=\\SI{{{row["C"]}}}{{ml/cmH_{2}O}}$, $R=\\SI{{{row["R"]}}}{{cmH_{2}O/l/s}}$, $p=\\SI{{{row["PINSP"]}}}{{cmH_{2}O}}$, r=$\\SI{{{row["RR"]}}}{{bpm}}$,$I:E=1:2$,$O_2=90-100$,$PEEP=\\SI{{{row["PEEP"]}}}{{cmH_{2}O}}$}}\n')
    toprint.append(f'%%%C={row["C"]}, R={row["R"]}, p={row["PINSP"]}, r={row["RR"]}, I:E=1:2, $O_2$=90-100, PEEP={row["PEEP"]}}}\n')
    toprint.append(f'C={row["C"]}, R={row["R"]}, p={row["PINSP"]}, r={row["RR"]}, PEEP={row["PEEP"]}}}')
    toprint.append(r'''
      \end{subfigure}'''
    )
  toprint.append(r'''  \caption{I:E=1:2, $O_2$=90-100}
    \end{figure}
  ''')

  toprint.append(r'''
    \end{document}
  ''')

  print(''.join(toprint))

if __name__ == '__main__':
  import argparse
  import matplotlib
  import style

  parser = argparse.ArgumentParser(description='present plots in LaTeX')
  parser.add_argument('input', help='plots to include', nargs='+')
  parser.add_argument('-l', '--logbook', type=str, help='logbook location', default='../Data/logbook.csv')
  parser.add_argument('-m', '--mhra', type=str, help='MHRA logbook location', default='../Data/logbook.MHRA.csv')
  args = parser.parse_args()

  regexp = re.compile('.*/(?P<datafile>\w+\.txt)_(?P<plotset>\d+)\..*') # accept PNG and PDF
  meta = read_meta_csv(args.logbook)
  mhra = read_mhra_csv(args.mhra)

  resistances = [5,20,10]
  compliances = [50,20,50]

  mhra = associate_plots(meta=meta, mhra=mhra, files=args.input, regexp=regexp, resistances=resistances, compliances=compliances)

  print_latex(mhra, 'Plots from Apr 3, 2020')
