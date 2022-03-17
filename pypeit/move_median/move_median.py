import warnings

try:
    from bottleneck import move_median as move_median_
    def move_median(seq, window):
        return move_median_(seq, window=window)
    #warnings.warn("using bottleneck")
except:
    warnings.warn('Unable to load bottleneck moving median.  Try reinstalling bottleneck.  In the '
                  'meantime, falling back on the slower pure python code.')
    from pypeit.move_median.mmpy import move_median
