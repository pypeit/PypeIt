from collections import deque
from bisect import insort, bisect_left
import itertools

def move_median(seq_pad, window_size):
    """
    Code contributed by Peter Otten, made to be consistent with
    scipy.ndimage.filters.median_filter by Joe Hennawi.

    See discussion at:
    http://groups.google.com/group/comp.lang.python/browse_thread/thread/d0e011c87174c2d0
    """
    odd = (window_size % 2) == 1

    seq_pad = iter(seq_pad)
    d = deque()
    s = []
    result = []
    for item in itertools.islice(seq_pad, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2] if odd else (s[len(d)//2] + s[len(d)//2 - 1])/2)
    m = window_size // 2
    for item in seq_pad:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m] if odd else (s[m] + s[m-1])/2)
    return result
