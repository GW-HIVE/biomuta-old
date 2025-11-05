from cmath import nan
import re

def convert_NA(NA_value):
    if NA_value in ['NA', 'None']:
        NA_value = nan
    return NA_value


def is_valid_enst(enst):
    # Check if the ENST matches the standard pattern and length
    return bool(re.match(r"^ENST\d{11}(\.\d+)?$", enst))