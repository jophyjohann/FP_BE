import numpy as np
import copy

def read_file(file_path):
    """Read file and return a directory with time, content, channels and counts"""
    input_file = open(file_path, 'r')
    # get the time
    time = int(input_file.readline()[-8:])
    # get the content
    content = np.loadtxt(file_path, skiprows=5, delimiter=',', usecols=(0, 1))

    # store time, content, channels, counts and statistical uncertainty of counts in a directory
    dataSet = {
        'name': file_path.split('/')[-1].split('.')[0],
        'time': time, 'content': content,
        'channel': np.transpose(content)[0],
        'counts': np.transpose(content)[1],
        'counts_uncert': np.maximum(np.sqrt(np.transpose(content)[1]), np.ones(len(np.transpose(content)[1])))
    }
    return dataSet


def subtract_file(dataSet1, dataSet2):
    """Subtract the counts of two data sets and calculates the resulting uncertainties.
    The different durations of measuring are considered.
    The resulting data set will have the 'time' of dataSet1"""
    # create a copy of dataSet1
    dataSet = copy.copy(dataSet1)
    dataSet_diff = []
    dataSet_uncert = []
    # scale the counts in dataSet2 according to the different durations of measuring
    scale_factor = float(dataSet1['time'])/float(dataSet2['time'])
    for i in range(0, len(dataSet1['counts'])):
        # minimum is 0 events and uncertainty of 1
        dataSet_diff.append(max(dataSet1['counts'][i] - dataSet2['counts'][i]*scale_factor, 0))
        dataSet_uncert.append(max(np.sqrt(dataSet1['counts_uncert'][i] ** 2 + (dataSet2['counts_uncert'][i] ** 2) * (scale_factor ** 2)), 1))
    dataSet['counts'] = np.array(dataSet_diff)
    dataSet['counts_uncert'] = np.array(dataSet_uncert)
    return dataSet


def rebin_file(dataSet1, rebin=1):
    """Rebins channel, counts, energy, etc. of a dataset. The number of merged bins is given by rebin"""
    # create a copy of dataSet1
    dataSet = copy.copy(dataSet1)
    # loop over all histograms
    for measurement in dataSet:
        # skip 'name', 'time', 'content'
        if measurement in ['name', 'time', 'content']:
            continue
        rebinned_array = []
        # loop over bins
        for i in range(0, int(len(dataSet[measurement]) / rebin) + 1):
            if i < int(len(dataSet[measurement]) / rebin):
                rebin_num = rebin
            else:
                # the remaining bins
                rebin_num = len(dataSet[measurement]) - i * rebin
                if rebin_num == 0:
                    continue
            if measurement in ['channel', 'energy']:
                rebinned_array.append(dataSet[measurement][i * rebin])
                continue
            value = 0
            for j in range(i * rebin, i * rebin + rebin_num):
                if not 'uncert' in measurement:
                    value +=  dataSet[measurement][j]
                else:
                    value = np.sqrt(value ** 2 + dataSet[measurement][j] ** 2)
            rebinned_array.append(value)

        rebinned_array = np.array(rebinned_array)
        dataSet[measurement] = rebinned_array
    return dataSet


def calibrate_dataSets(dataSets, cal_parameters):
    for dataSet in dataSets:
        dataSet['energy'] = cal_parameters[0] * dataSet['channel'] + cal_parameters[1]