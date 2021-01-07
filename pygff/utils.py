import sys
import gzip
import time


class ProgressBar(object):

    def __init__(self, total, prefix="", suffix="", ncol=60, file=sys.stderr):
        self.count = total
        self.prefix = prefix
        self.suffix = suffix
        self.ncol = ncol
        self.file = file
        self.done = 0
        self.last_time = time.time()

    def update(self, amount):
        self.done += amount
        x = int(self.ncol*self.done/self.count)
        self.file.write("%s[%s%s] %i/%i %s\r" % (
            self.prefix, "#"*x, " "*(self.ncol-x), self.done, self.count,
            self.suffix))
        self.file.flush()

    def __enter__(self):
        self.update(0)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.write("%s[%s%s] %i/%i %s\r\n" % (
            self.prefix, "#"*self.ncol, " "*0, self.done, self.count,
            self.suffix))
        self.file.flush()
