def join_unique(series):
    """
    This function takes a pandas Series, filters out None values, and returns a string of unique values joined by a comma.

    Parameters:
    - series: A pandas Series object.

    Returns:
    - A string of unique values in the series, joined by a comma. If a value in the series is None, it is not included in the output string.
    """
    # Filter out None values before converting to set and joining
    filtered_series = [str(item) for item in series if item is not None]
    unique_items = set(filtered_series)
    return ', '.join(unique_items)


