from scipy.interpolate import interp1d
import numpy as np

def series_interp(data = []):
    y = [k for k in data if k >= 0]
    x = [i for i in range(len(data)) if data[i] >= 0 ]
    if 0 not in x:
        y = [0] + y
        x = [0] + x
    f = interp1d(x,y,kind='quadratic')
    return [data[i] if data[i] >= 0 else f(i) for i in range(len(data))]

def window_avg(data = [], window = 7):
    avg = []
    for i in range(len(data)):
        if i < window:
            avg.append(np.mean(data[:i+1]))
        else:
            avg.append(np.mean(data[i-window+1:i+1]))
    return avg

def run_avg(x, wind=1):
    avg = []
    for i in range(len(x)):
        if i >= np.floor(wind/2) and i < np.ceil(len(x) - wind/2):
            avg.append(np.mean(x[int(i+np.floor(-wind/2+1)):int(i+np.floor(wind/2+1))]))
        elif i < wind/2:
            avg.append(np.mean(x[:int(i+wind/2+1)]))
        else:
            avg.append(np.mean(x[int(i-wind/2+1):]))
    return avg


def make_fig_dates():
    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)
    formatter = mdates.DateFormatter("%Y-%m-%d")
    locator = mdates.MonthLocator()
    locator_days = mdates.DayLocator()
    ax.xaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_locator(locator)
    ax.set_xlabel("Date",fontsize=20)
    ax.xaxis.set_ticks(tick_dates)
    ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylabel("Percent of Population (%)",fontsize=20)
    ax.tick_params(labelsize=15, direction='out', length=6, width=2)
    #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
    ax.yaxis.offsetText.set_fontsize(20)
    ax.set_prop_cycle(color=[
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    '#17becf', '#9edae5'])
    
    return f,ax

def make_fig(f=None,ax=None):
    if f == None and ax == None:
        f = plt.figure(figsize=(10,10))
        ax = f.add_subplot(111)
    ax.set_xlabel("Frequency",fontsize=18)
    ax.set_ylabel("Percent of Population",fontsize=18)
    def func(x,pos):
        if str(np.around(x,2)) ==str(x):
            return x
        return np.format_float_scientific(x,precision=2)
    ax.tick_params(labelsize=16, direction='out', length=6, width=2)
    # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)
    ax.yaxis.offsetText.set_fontsize(14)
    ax.set_prop_cycle(color=[
    '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
    '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
    '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
    '#17becf', '#9edae5'])
    return f,ax


def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])
