
DATADIR='/Users/sviel/MVM/Data'
#python python/combine.py ${DATADIR} -p --db-google-id 1dC59lsoh9G72JBDVbFGRUV0WQDa38KB_   # Old spreadsheet from Shawn, instead use new default
#python python/combine.py ${DATADIR} -p --db-range-name "20200408 ISO!A2:AP4"              # Success, plots the first metadata row
#python python/combine.py ${DATADIR} -p --db-range-name "20200412 ISO!A2:AP15" --campaign "Run_15"

## working command from Paolo
#python python/combine.py ${DATADIR} -p --mvm-col='mvm_col_arduino' -d plots_iso_Apr13c -f VENTILATOR_12042020_CONTROLLED_FR12_PEEP10_PINSP35_C20_R50_RATIO025.txt

## ignore_sim doesn't work
#python python/combine.py ${DATADIR} -p --mvm-col='mvm_col_arduino' -d plots_iso_Apr13_ignoresim --ignore_sim -f VENTILATOR_12042020_CONTROLLED_FR20_PEEP5_PINSP30_C20_R20_RATIO050_leak.txt

## plot TRIUMF data
#python python/triumf.py ${DATADIR} -p -d plots_triumf_Apr10
python python/triumf.py ${DATADIR} -p -d plots_triumf_Apr10b --db-google-id 1aQjGTREc9e7ScwrTQEqHD2gmRy9LhDiVatWznZJdlqM



#python latexify.py ../Plots/Run\ 9\ Apr\ 3\ 2020/*txt*pdf > ../Plots/Run\ 9\ Apr\ 3\ 2020/summary.tex
#python get_tables.py plots_iso/*json --output-dir=plots_iso
