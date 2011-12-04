## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

#############################################################################
# 
# Histogram class: 
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr> 
# last revision: 1999-7-6
#
# HistogramRI class written by Ruth Huey
# HistogramRI adds a reverse Index to Histogram class 
#
#############################################################################


import numpy.oldnumeric as Numeric; N = Numeric

class Histogram:

    """Histogram in one variable

    Constructor: Histogram(|data|, |bins|, |range|=None)

    Arguments:

    |data| -- a sequence of data points

    |bins| -- the number of bins into which the data is to be sorted

    |range| -- a tuple of two values, specifying the lower and
               the upper end of the interval spanned by the bins.
               Any data point outside this interval will be ignored.
               If no range is given, the smallest and largest
               data values are used to define the interval.

    The number of points in a bin can be obtained by indexing the
    histogram with the bin number. Application of len() yields the
    number of bins. A histogram thus behaves like a sequence of
    numbers.
    """

    def __init__(self, data, nbins, range=None):
        if range is None:
            self.min = N.minimum.reduce(data)
            self.max = N.maximum.reduce(data)
        else:
            self.min, self.max = range
        self.min = self.min+0.
        self.max = self.max+0.
        self.bin_width = (self.max-self.min)/nbins
        if self.bin_width==0:
            print 'range is 0 so set bin_width to 1.'
            self.bin_width = 1.
        self.array = N.zeros((nbins, 2), N.Float)
        self.array[:, 0] = self.min + self.bin_width*(N.arange(nbins)+0.5)
        self.addData(data)

    def __len__(self):
        return self.array.shape[0]

    def __getitem__(self, index):
        return self.array[index]

    def __getslice__(self, first, last):
        return self.array[first:last]

    def addData(self, data):
        """Add the values in |data| (a sequence of numbers) to the
        originally supplied data. Note that this does not affect the
        default range of the histogram, which is fixed when the
        histogram is created.
        """
        n = (len(data)+999)/1000
        for i in range(n):
            self._addData(data[1000*i:1000*(i+1)])

    def _addData(self, data):
        data = N.array(data, N.Float)
        data = N.repeat(data, N.logical_and(N.less_equal(data, self.max),
                                            N.greater_equal(data, self.min)))
        data = N.floor((data - self.min)/self.bin_width).astype(N.Int)
        self.rIdata = data
        nbins = self.array.shape[0]
        histo = N.add.reduce(N.equal(N.arange(nbins)[:,N.NewAxis], data), -1)
        histo[-1] = histo[-1] + N.add.reduce(N.equal(nbins, data))
        self.array[:, 1] =  self.array[:, 1] + histo

    def normalize(self, norm=1.):
        "Scales all counts by the same factor such that their sum is |norm|."
        self.array[:, 1] = norm*self.array[:, 1]/N.add.reduce(self.array[:, 1])


class HistogramRI(Histogram):
    """ This class is based on Histogram class developed by K.Hinsen. It adds
    a method, 'createReverseIndex', which builds a list of data points in each
    bin"""

    def __init__(self, data, nbins, range=None):
        Histogram.__init__(self, data, nbins, range)
        self.data = data

    
    def createReverseIndex(self):
        #first build the data in the correct form:
        N = Numeric
        nbins = self.array.shape[0]
        data = N.array(self.data, N.Float)
        data = N.repeat(data, N.logical_and(N.less_equal(data, self.max),
                    N.greater_equal(data, self.min)))
        data = N.floor((data-self.min)/self.bin_width).astype(N.Int)
        self.cRI_data = data
        #check for pt outside of nbins
        #sometimes max point falls outside of right-most bin
        #first stuff works with python2.0 second with both
        try:
            self.cRI_data = N.putmask(data, N.greater_equal(data,nbins),nbins-1)
        except:
            for i in range(data.shape[0]):
                if data[i]>=nbins: data[i] = nbins-1
        self.reverseIndex = []
        for i in range(nbins):
            newentry = N.nonzero(N.equal(data, i)).tolist()
            #print i,'th entry =', newentry
            self.reverseIndex.append(newentry)
            

if __name__ == '__main__':


    data = N.arange(500)
    #data = N.arange(5000)
    #h = Histogram(data, 10)
    hRI = HistogramRI(data, 10)
