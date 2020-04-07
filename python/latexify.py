''' Put plots together in LaTeX documents '''

import re
from db import *

def associate_plots(meta, files, regexp, table_entry='MVM_filename'):
  # loop over available files
  meta['fname'] = ''

  for fname in files:
    m = regexp.match(fname)
    if m:
      # TODO: we are ignoring the plotset info!
      plotset = int(m.group('plotset'))
      if plotset != 0:
        raise RuntimeError('Did plot naming conventions change? check and fix the code')
      datafile = f'{m.group("datafile")}.txt'
      sel = meta[table_entry]==datafile # we select the row corresponding to this datafile

      meta.loc[sel, 'fname'] = fname
      meta.loc[sel, 'plot'] = 1
    else:
      raise RuntimeError(f'Unable to parse the name of the MVM file {fname}')

  return meta

def print_latex(meta, title, figure_caption="", rows=3, columns=2):
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

  for i, row in meta.iterrows():
    newfig = i % (rows * columns) == 0
    fname_form = row["fname"].split('/')[-1] # relative path to current dir

    if newfig:
      if i != 0:
        toprint.append(r'''  \caption{''' + figure_caption + r'''}
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
    toprint.append(f'%%% $TV=\\SI{{{row["TV"]}}}{{ml}}$, $C=\\SI{{{row["C"]}}}{{ml/cmH_{2}O}}$, $R=\\SI{{{row["R"]}}}{{cmH_{2}O/l/s}}$, $p=\\SI{{{row["plateau"]}}}{{cmH_{2}O}}$, r=$\\SI{{{row["rate"]}}}{{bpm}}$,$I:E={row["ratio"]}$, $O_2={row["O2"]}$, $PEEP=\\SI{{{row["PEEP"]}}}{{cmH_{2}O}}$}}\n')
    toprint.append(f'%%%TV={row["TV"]}, C={row["C"]}, R={row["R"]}, p={row["plateau"]}, r={row["rate"]}, I:E={row["ratio"]}, $O_2={row["O2"]}$, PEEP={row["PEEP"]}}}\n')
    toprint.append(f'TV={row["TV"]}, C={row["C"]}, R={row["R"]}, p={row["plateau"]}, r={row["rate"]}, PEEP={row["PEEP"]}}}, I:E={row["ratio"]}')
    toprint.append(r'''
      \end{subfigure}'''
    )
  toprint.append(r'''  \caption{''' + figure_caption + r'''}
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
  parser.add_argument("--db-google-id", type=str, help="name of the Google spreadsheet ID for metadata", default="1aQjGTREc9e7ScwrTQEqHD2gmRy9LhDiVatWznZJdlqM")
  parser.add_argument("--db-range-name", type=str, help="name of the Google spreadsheet range for metadata", default="20200407 ISO!A2:AZ")
  args = parser.parse_args()

  test_data = read_online_spreadsheet(spreadsheet_id=args.db_google_id, range_name=args.db_range_name)

  regexp = re.compile('.*/(?P<datafile>.+)_(?P<plotset>\d+)\.\w\w\w') # accept PNG and PDF

  test_data = associate_plots(meta=test_data, files=args.input, regexp=regexp, table_entry='MVM_filename')

  print_latex(test_data, title='Plots from Apr 7, 2020', figure_caption='$O_2=21\\%$')
