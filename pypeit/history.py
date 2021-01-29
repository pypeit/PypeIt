
from astropy.time import Time

class History:
    def __init__(self, header=None):
        self.history = []
        if header is not None and 'HISTORY' in header:
            for card in header['HISTORY']:
                self.history.append(str(card))

    def append(self, history, add_date=True):
        if add_date:
            self.history.append(f'{Time.now().to_value("isot", subfmt="date_hm")} {history}')
        else:
            self.history.append(history)
        
    def write_to_header(self, header):
        for line in self.history:
            header['history'] = line

