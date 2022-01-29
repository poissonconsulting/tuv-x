#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

from yahoo_subs import load_csv_data
from plot_subs import line_plot, scat_plot, line_scat_plot
from pltClass import FinInst
import numpy as np
import sys
from datetime import datetime as dt
import matplotlib.pyplot as plt

"""
The main program
"""

#--------------------------
# get start,end dates
#--------------------------
dat_prompt = '\nEnter start,end dates(Jun 1,1960 -> Jun 3,1961): '

try:
   dates = input(dat_prompt)
except EOFError:
   print('start, end dates not input')
   sys.exit(-1)

dates = dates.strip()
datem = dates.replace(' ','')
datem = datem.split('->')
dates = dates.split('->')

sdate = dt.strptime(datem[0],"%b%d,%Y")
edate = dt.strptime(datem[1],"%b%d,%Y")
print(f"Start date obj = {sdate}")
print(f"End   date obj = {edate}")

if( edate < sdate ):
   print(f"End date {dates[1]} is before start date {dates[0]}")
   sys.exit(-1)

dates = [sdate.strftime("%d-%m-%Y"), edate.strftime("%d-%m-%Y")]

plt.style.use('seaborn')

finInst = FinInst()

while True:
#--------------------------
# get plot parameters
#--------------------------
   if( finInst.getPltParms() == True ): break
  
   symbol = finInst.dict['sym']
   try:
      findata = load_csv_data(symbol,day_begin=dates[0],day_end=dates[1])
   except IndexError:
       print(f"\nCould not acquire data for {symbol}")
       continue

   if( len(findata) < 2 ):
       print(f"\nCould not acquire data for {symbol}")
       continue

   isIndex = True if symbol[0] == '^' else False

   seriesDate = []
   seriesClose = []
   seriesVol = []

   Skip = False
   for line in findata[1:]:
      tokens = line.split(',')
      if( 'null' in tokens ):
         print(f"Null entries in {symbol} data")
         Skip = True
         break
#--------------------------
# extract date,
# adjusted closing price,
# and volume
#--------------------------
      seriesDate.append(tokens[0])
      seriesClose.append(float(tokens[-2]))
      seriesVol.append(float(tokens[-1]))

   if( Skip ): continue

   npDate = [dt.strptime(date,'%Y-%m-%d').date() for date in seriesDate]

   npClose = np.array(seriesClose)
   npVol   = np.array(seriesVol)

#--------------------------
# price or volume data?
#--------------------------
   if( finInst.dict['inp'] == 'price' ):
      Data = npClose
      Title = symbol + ' Close'
   elif( finInst.dict['inp'] == 'volume' ):
      Data = npVol
      Title = symbol + ' Volume'

#--------------------------
# which function applied to data?
#--------------------------
   function = finInst.dict['ftn']
   dateSndx = 1 if function == 'reldel' else 0
   if( function == 'stkdel' ):
      Series = 100.*(Data[:] / Data[0] - 1.)
      yLabel = '% delta'
      yLabel += ' Index' if isIndex else ' Price'
   elif( function == 'none' ):
      Series = Data
      yLabel = 'Closing'
      yLabel += ' Index' if isIndex else ' Price'
   elif( function == 'reldel' ):
      dataLndx = len(Data)
      Series = 100.*(Data[1:dataLndx] / Data[0:dataLndx-1] - 1.)
      yLabel = '% delta'
      yLabel += ' Index' if isIndex else ' Price'

   Dates = npDate[dateSndx:]
#--------------------------
# now make a simple plot
#--------------------------
   if( finInst.dict['plt'] == 'line' ):
      line_plot(symbol,yLabel,Title,Dates,Series)
   elif( finInst.dict['plt'] == 'scatter' ):
      scat_plot(symbol,yLabel,Title,Dates,Series)
   elif( finInst.dict['plt'] == 'line&scatter' ):
      line_scat_plot(symbol,yLabel,Title,Dates,Series)

print(f"\n")
