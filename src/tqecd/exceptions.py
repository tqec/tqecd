# N818 is not applied here: renaming this public exception would break users.
class TQECDException(Exception):  # noqa: N818
    pass


class TQECDWarning(Warning):
    pass
