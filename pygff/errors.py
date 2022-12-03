class SeqExtractError(Exception):
    pass

class PositionNotSpecified(SeqExtractError):
    pass

class ChromosomeNotSpecified(SeqExtractError):
    pass