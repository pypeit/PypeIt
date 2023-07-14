import warnings

try:
    from bottleneck import move_median as move_median_
    def move_median(seq, window):
        # TODO JFH This attempts to make the behavior of bottleneck moving_median and python moving median equivalent.
        # Perhaps we should be crashing here instead of returning the input? fast_running_median issues a warning but
        # does not fault
        return seq if len(seq) == 0 else move_median_(seq, window=window)
    #warnings.warn("using bottleneck")
except:
    warnings.warn('Unable to load bottleneck moving median.  Try reinstalling bottleneck.  In the '
                  'meantime, falling back on the slower pure python code.')
    from pypeit.move_median.mmpy import move_median
