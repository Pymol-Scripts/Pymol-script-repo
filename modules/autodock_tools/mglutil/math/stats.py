import warnings

def stats(values):
    """returns the mimn, max, mean and standard deviation of a list of values"""

    warnings.warn(
        "\n\nWARNING!! This function has been deprecated!!\n \
        Use the stats in Volume/Grid3D.py\n",DeprecationWarning,2)
    
    npts = len(values)
    if npts:
        from math import sqrt
        sum = 0.0
        sumsq = 0.0
        mini = maxi = values[0]
        for v in values:
            sum += v
            sumsq += float(v)*float(v)
            if v<mini:
                mini = v
            if v>maxi:
                maxi = v
        mean = float(sum)/npts
        stdev = sqrt(( sumsq - (sum*sum/float(npts)))/(npts-1))
        return mini, maxi, mean, stdev

    else:
      return (0., 0., 1., 1.)
   
