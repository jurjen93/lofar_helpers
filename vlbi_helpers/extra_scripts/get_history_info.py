import os
import re

def isfloat(num):
    """
    Check if value is a float
    """
    try:
        float(num)
        return True
    except ValueError:
        return False

def parse_history(ms, hist_item):
    """
    Grep specific history item from MS

    :param ms: measurement set
    :param hist_item: history item

    :return: parsed string
    """
    hist = os.popen('taql "SELECT * FROM '+ms+'::HISTORY" | grep '+hist_item).read().split(' ')
    for item in hist:
        if hist_item in item and len(hist_item)<=len(item):
            return item
    print('WARNING:' + hist_item + ' not found')
    return None


def get_time_preavg_factor(ms):
    """
    Get pre-time averaging factor (given by demixer.timestep)

    :param ms: measurement set

    :return: averaging integer
    """
    parse_str = "demixer.timestep="
    parsed_history = parse_history(ms, parse_str)
    avg_num = re.findall(r'\d+', parsed_history.replace(parse_str, ''))[0]
    if avg_num.isdigit():
        return int(float(avg_num))
    elif isfloat(avg_num):
        print("WARNING: parsed factor is not a digit but a float")
        return float(avg_num)
    else:
        print("WARNING: parsed factor is not a float or digit")
        return None